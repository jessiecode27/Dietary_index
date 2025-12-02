#' ACS2020_NHANES_FPED
#'
#' Calculate the American Cancer Society 2020 dietary index (ACS2020) for the NHANES_FPED data (after 2005) within 1 step
#' @import dplyr
#' @import readr
#' @import haven
#' @param FPED_IND_PATH The file path for the FPED IND data for day 1. The file name should be like: fped_dr1iff.sas7bdat
#' @param NUTRIENT_IND_PATH The file path for the NUTRIENT IND data for day 1. The file name should be like: DR1IFF_J
#' @param FPED_IND_PATH2 The file path for the FPED IND data for day 2. The file name should be like: fped_dr2iff.sas7bdat
#' @param NUTRIENT_IND_PATH2 The file path for the NUTRIENT IND data for day 2. The file name should be like: DR2IFF_J
#' @param SSB_code The food code for sugar sweetened beverage, default is the SSB code from 17-18 FNDDS file.
#' @return The ACS2020 and its component scores and serving sizes
#' @examples
#' data("NHANES_20172018")
#' ACS2020_NHANES_FPED(FPED_IND_PATH = NHANES_20172018$FPED_IND, NUTRIENT_IND_PATH = NHANES_20172018$NUTRIENT_IND, FPED_IND_PATH2 = NHANES_20172018$FPED_IND2, NUTRIENT_IND_PATH2 = NHANES_20172018$NUTRIENT_IND2)
#' @export

ACS2020_NHANES_FPED = function(FPED_IND_PATH = NULL, NUTRIENT_IND_PATH = NULL, FPED_IND_PATH2 = NULL, NUTRIENT_IND_PATH2 = NULL, SSB_code = NULL) {
    # stop if the input data is not provided for any day
    if (is.null(FPED_IND_PATH) & is.null(NUTRIENT_IND_PATH) & is.null(FPED_IND_PATH2) & is.null(NUTRIENT_IND_PATH2)) {
        stop("Please provide the file path for the FPED and NUTRIENT data, day 1 or day 2 or day 1 and day 2.")
    }

    if (is.null(SSB_code)) {
        # load the SSB codes from 17-18 FNDDS file as default
        data("SSB_FNDDS_1718")
        SSB = unique(SSB_FNDDS_1718$`Food code`)
        message("Since no SSB code is provided, the default SSB code from 17-18 FNDDS file is used.")
    } else {
        # use the provided SSB code
        SSB = SSB_code
    }

    # if only day 1 data is provided
    if (!is.null(FPED_IND_PATH) & !is.null(NUTRIENT_IND_PATH)) {
        if (is.character(FPED_IND_PATH) == TRUE) {
            FPED_IND = read_sas(FPED_IND_PATH)
        } else {
            FPED_IND = FPED_IND_PATH
        }

        if (is.character(NUTRIENT_IND_PATH) == TRUE) {
            NUTRIENT_IND = read_xpt(NUTRIENT_IND_PATH)
        } else {
            NUTRIENT_IND = NUTRIENT_IND_PATH
        }

        if (!("DR1ILINE" %in% colnames(FPED_IND)) | !("DR1ILINE" %in% colnames(NUTRIENT_IND))) {
            stop("Please use individual-level first day data. Individual-level nutrient data should be like DR1IFF_J.XPT. Individual-level FPED data should be like fped_dr1iff_1718.sas7bdat")
        }

        # Select only the high quality data
        NUTRIENT_IND = NUTRIENT_IND %>%
            filter(DR1DRSTZ == 1) %>%
            arrange(SEQN)

        FPED_IND = FPED_IND %>%
            arrange(SEQN)

        # merge the FPED and nutrient data. Here, the data is all the individual consumed foods and nutrients for each individual (SEQN). The DR1IFDCD.x is the food code for the consumed food.
        COHORT = NUTRIENT_IND %>%
            left_join(FPED_IND, by = c("SEQN", "DR1ILINE"))

        # Create the serving size variables for ACS2020 calculation
        COHORT_summary = COHORT %>%
            dplyr::mutate(
                # create the variable for added sugars from SSB
                ADDED_SUGAR_SSB_SERV = case_when(
                    DR1IFDCD.x %in% SSB ~ DR1I_ADD_SUGARS,
                    TRUE ~ 0
                )
            ) %>%
            # group by individual
            dplyr::group_by(SEQN) %>%
            # summarize to get the total servings of each food group for each individual
            dplyr::summarize(
                ENERGY = sum(DR1IKCAL),
                RIAGENDR = min(RIAGENDR),
                VEG_SERV = sum(DR1I_V_DRKGR + (DR1I_V_REDOR_TOTAL + DR1I_V_OTHER + DR1I_V_STARCHY_OTHER) / 0.5),
                # count the unique number of DR1IFDCD.x with any values in vegetables (DR1I_V_DRKGR, DR1I_V_REDOR_TOTAL, DR1I_V_OTHER, DR1I_V_STARCHY_OTHER)
                VEG_ITEMS_SERV = length(unique(DR1IFDCD.x[DR1I_V_DRKGR > 0 | DR1I_V_REDOR_TOTAL > 0 | DR1I_V_OTHER > 0 | DR1I_V_STARCHY_OTHER > 0])),
                FRT_SERV = sum(DR1I_F_TOTAL - DR1I_F_JUICE),
                # count the unique number of DR1IFDCD.x with any values in fruits but not fruit juice (DR1I_F_TOTAL - DR1I_F_JUICE)
                FRT_ITEMS_SERV = length(unique(DR1IFDCD.x[DR1I_F_TOTAL > 0 & DR1I_F_JUICE == 0])),
                WGRAIN_SERV = sum(DR1I_G_WHOLE / 0.035274),
                REDPROC_MEAT_SERV = sum((DR1I_PF_CUREDMEAT / 1.5) + ((DR1I_PF_MEAT + DR1I_PF_ORGAN) / 4)),
                # The daily servings of highly processed foods and refined grains per 1000 kcal
                # the ultra-processed variable for the score should not double count foods included in other parts of the score, for example, sugar-sweetened beverages or processed meats
                # In NHANES, this should be refined grains
                HPFRG_SERV = sum(DR1I_G_REFINED / 0.035274),
                SSB_FRTJ_SERV = sum((ADDED_SUGAR_SSB_SERV * 4 / 26))
            )

        # select only participants with more than 0 kcal intake
        COHORT_summary = COHORT_summary %>%
            filter(ENERGY > 0)

        ## ACS2020 score calculation
        COHORT_summary_ACS2020 = ACS2020_V2(
            SERV_DATA = COHORT_summary,
            RESPONDENTID = COHORT_summary$SEQN,
            GENDER = COHORT_summary$RIAGENDR,
            TOTALKCAL_ACS2020 = COHORT_summary$ENERGY,
            VEG_SERV_ACS2020 = COHORT_summary$VEG_SERV,
            VEG_ITEMS_SERV_ACS2020 = COHORT_summary$VEG_ITEMS_SERV,
            FRT_SERV_ACS2020 = COHORT_summary$FRT_SERV,
            FRT_ITEMS_SERV_ACS2020 = COHORT_summary$FRT_ITEMS_SERV,
            WGRAIN_SERV_ACS2020 = COHORT_summary$WGRAIN_SERV,
            REDPROC_MEAT_SERV_ACS2020 = COHORT_summary$REDPROC_MEAT_SERV,
            HPFRG_SERV_ACS2020 = COHORT_summary$HPFRG_SERV,
            SSB_FRTJ_SERV_ACS2020 = COHORT_summary$SSB_FRTJ_SERV
        )

        # rename RESPONDENTID to SEQN and GENDER to RIAGENDR
        COHORT_summary_ACS2020 = COHORT_summary_ACS2020 %>%
            rename(SEQN = RESPONDENTID, RIAGENDR = GENDER)

    }

    # if only day 2 data is provided
    if (!is.null(FPED_IND_PATH2) & !is.null(NUTRIENT_IND_PATH2)) {
        if (is.character(FPED_IND_PATH2) == TRUE) {
            FPED_IND2 = read_sas(FPED_IND_PATH2)
        } else {
            FPED_IND2 = FPED_IND_PATH2
        }

        if (is.character(NUTRIENT_IND_PATH2) == TRUE) {
            NUTRIENT_IND2 = read_xpt(NUTRIENT_IND_PATH2)
        } else {
            NUTRIENT_IND2 = NUTRIENT_IND_PATH2
        }

        if (!("DR2ILINE" %in% colnames(FPED_IND2)) | !("DR2ILINE" %in% colnames(NUTRIENT_IND2))) {
            stop("Please use individual-level second day data. Individual-level nutrient data should be like DR2IFF_J.XPT. Individual-level FPED data should be like fped_dr2iff_1718.sas7bdat")
        }


        NUTRIENT_IND2 = NUTRIENT_IND2 %>%
            filter(DR2DRSTZ == 1) %>%
            arrange(SEQN)

        FPED_IND2 = FPED_IND2 %>%
            arrange(SEQN)

        COHORT2 = NUTRIENT_IND2 %>%
            left_join(FPED_IND2, by = c("SEQN", "DR2ILINE"))

        COHORT2_summary = COHORT2 %>%
            dplyr::mutate(
                # create the variable for added sugars from SSB
                ADDED_SUGAR_SSB_SERV = case_when(
                    DR2IFDCD.x %in% SSB ~ DR2I_ADD_SUGARS,
                    TRUE ~ 0
                )
            ) %>%
            # group by individual
            dplyr::group_by(SEQN) %>%
            # summarize to get the total servings of each food group for each individual
            dplyr::summarize(
                ENERGY = sum(DR2IKCAL),
                RIAGENDR = min(RIAGENDR),
                VEG_SERV = sum(DR2I_V_DRKGR + (DR2I_V_REDOR_TOTAL + DR2I_V_OTHER + DR2I_V_STARCHY_OTHER) / 0.5),
                VEG_ITEMS_SERV = length(unique(DR2IFDCD.x[DR2I_V_DRKGR > 0 | DR2I_V_REDOR_TOTAL > 0 | DR2I_V_OTHER > 0 | DR2I_V_STARCHY_OTHER > 0])),
                FRT_SERV = sum(DR2I_F_TOTAL - DR2I_F_JUICE),
                FRT_ITEMS_SERV = length(unique(DR2IFDCD.x[DR2I_F_TOTAL > 0 & DR2I_F_JUICE == 0])),
                WGRAIN_SERV = sum(DR2I_G_WHOLE / 0.035274),
                REDPROC_MEAT_SERV = sum((DR2I_PF_CUREDMEAT / 1.5) + ((DR2I_PF_MEAT + DR2I_PF_ORGAN) / 4)),
                # The daily servings of highly processed foods and refined grains per 1000 kcal
                # the ultra-processed variable for the score should not double count foods included in other parts of the score, for example, sugar-sweetened beverages or processed meats
                # In NHANES, this should be refined grains
                HPFRG_SERV = sum(DR2I_G_REFINED / 0.035274),
                SSB_FRTJ_SERV = sum((ADDED_SUGAR_SSB_SERV * 4 / 26))
            )

        # select only participants with more than 0 kcal intake
        COHORT2_summary = COHORT2_summary %>%
            filter(ENERGY > 0)

        ## ACS2020 calculation
        COHORT2_summary_ACS2020 = ACS2020_V2(
            SERV_DATA = COHORT2_summary,
            RESPONDENTID = COHORT2_summary$SEQN,
            GENDER = COHORT2_summary$RIAGENDR,
            TOTALKCAL_ACS2020 = COHORT2_summary$ENERGY,
            VEG_SERV_ACS2020 = COHORT2_summary$VEG_SERV,
            VEG_ITEMS_SERV_ACS2020 = COHORT2_summary$VEG_ITEMS_SERV,
            FRT_SERV_ACS2020 = COHORT2_summary$FRT_SERV,
            FRT_ITEMS_SERV_ACS2020 = COHORT2_summary$FRT_ITEMS_SERV,
            WGRAIN_SERV_ACS2020 = COHORT2_summary$WGRAIN_SERV,
            REDPROC_MEAT_SERV_ACS2020 = COHORT2_summary$REDPROC_MEAT_SERV,
            HPFRG_SERV_ACS2020 = COHORT2_summary$HPFRG_SERV,
            SSB_FRTJ_SERV_ACS2020 = COHORT2_summary$SSB_FRTJ_SERV
        )

        # rename RESPONDENTID to SEQN and GENDER to RIAGENDR
        COHORT2_summary_ACS2020 = COHORT2_summary_ACS2020 %>%
            rename(SEQN = RESPONDENTID, RIAGENDR = GENDER)
    }

    if (!is.null(FPED_IND_PATH) & !is.null(NUTRIENT_IND_PATH) & is.null(FPED_IND_PATH2) & is.null(NUTRIENT_IND_PATH2)) {
        message("Trans fat is not avaiable for NHANES, so it is not included in the ACS2020 score.")
        return(COHORT_summary_ACS2020)
    }

    if (is.null(FPED_IND_PATH) & is.null(NUTRIENT_IND_PATH) & !is.null(FPED_IND_PATH2) & !is.null(NUTRIENT_IND_PATH2)) {
        message("Trans fat is not avaiable for NHANES, so it is not included in the ACS2020 score.")
        return(COHORT2_summary_ACS2020)
    }

    # merge two days data if they both exist by creating columns with the same name but taking average of the original column
    if (!is.null(FPED_IND_PATH) & !is.null(NUTRIENT_IND_PATH) & !is.null(FPED_IND_PATH2) & !is.null(NUTRIENT_IND_PATH2)) {
        COHORT12_summary_ACS2020 = inner_join(COHORT_summary_ACS2020, COHORT2_summary_ACS2020, by = "SEQN") %>%
            mutate(
                # ACS2020_V2_ALL TOTALKCAL_ACS2020 ACS2020_VEG ACS2020_VEG_ITEMS ACS2020_FRT ACS2020_FRT_ITEMS ACS2020_WGRAIN ACS2020_REDPROC_MEAT ACS2020_HPFRG ACS2020_SSB_FRTJ
                ACS2020_V2_ALL = (ACS2020_V2_ALL.x + ACS2020_V2_ALL.y) / 2,
                RIAGENDR = RIAGENDR.x,
                TOTALKCAL_ACS2020 = (TOTALKCAL_ACS2020.x + TOTALKCAL_ACS2020.y) / 2,
                ACS2020_VEG = (ACS2020_VEG.x + ACS2020_VEG.y) / 2,
                ACS2020_VEG_ITEMS = (ACS2020_VEG_ITEMS.x + ACS2020_VEG_ITEMS.y) / 2,
                ACS2020_FRT = (ACS2020_FRT.x + ACS2020_FRT.y) / 2,
                ACS2020_FRT_ITEMS = (ACS2020_FRT_ITEMS.x + ACS2020_FRT_ITEMS.y) / 2,
                ACS2020_WGRAIN = (ACS2020_WGRAIN.x + ACS2020_WGRAIN.y) / 2,
                ACS2020_REDPROC_MEAT = (ACS2020_REDPROC_MEAT.x + ACS2020_REDPROC_MEAT.y) / 2,
                ACS2020_HPFRG = (ACS2020_HPFRG.x + ACS2020_HPFRG.y) / 2,
                ACS2020_SSB_FRTJ = (ACS2020_SSB_FRTJ.x + ACS2020_SSB_FRTJ.y) / 2
            ) %>%
            select(
                SEQN, RIAGENDR, ACS2020_V2_ALL, TOTALKCAL_ACS2020, ACS2020_VEG, ACS2020_VEG_ITEMS, ACS2020_FRT, ACS2020_FRT_ITEMS, ACS2020_WGRAIN, ACS2020_REDPROC_MEAT, ACS2020_HPFRG, ACS2020_SSB_FRTJ
            )

        return(COHORT12_summary_ACS2020)
    }
}
