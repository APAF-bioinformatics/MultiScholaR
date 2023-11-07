if(!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(testthat)
p_load(tidyverse)
p_load(ProteomeRiver)


test_check("ProteomeRiver")

## Load the packages






## clean_isoform_number

context("clean_isoform_number")


test_that( "clean_isoform_number function ", {expect_that(cleanIsoformNumber("Q8K4R4-2"), equals("Q8K4R4")) })




##--------

context("get_max_prob_future_map")


# phosphopeptide <- "AAAAAAAAAAGDS(0.794)DS(0.206)WDADTFSMEDPVRK"
test_that( "get_max_prob_future_map function ",  {expect_that(getMaxProbFutureMap("AAAAAAAAAAGDS(0.794)DS(0.206)WDADTFSMEDPVRK", 1)[[1]], equals(0.794))} )
test_that( "get_max_prob_future_map function ",  {expect_that(getMaxProbFutureMap("IQAAAS(1)PPANATAASDTNAGDR", 1)[[1]], equals(1))} )
test_that( "get_max_prob_future_map function ",  {expect_that(getMaxProbFutureMap("IQAAAS(1)PPANAT(0.99)AASDTNAGDR", 2), equals(list(c(1, 0.99) ) ))} )
test_that( "get_max_prob_future_map function ",  {expect_that(getMaxProbFutureMap("AAS(1)PQALT(0.889)PT(0.101)LALT(0.01)LPP", 2)[[1]], equals(c(1.000, 0.889)  ))} )

peptides <- list( "AAAAT(0.003)APPS(0.997)PGPAQPGPR",
                  "AAAAT(0.013)APPS(0.987)PGPAQPGPR"	)
test_that( "get_max_prob_future_map function ",  {expect_that(getMaxProbFutureMap(peptides, list(1, 1)) ,
                                                              equals( list( 0.997, 0.987) ))} )





##--------

context("get_best_position")

test_that( "get_best_position function ",  {expect_that(getBestPosition("AAAAAAAAAAGDS(0.794)DS(0.206)WDADTFSMEDPVRK") ,
                                                        equals( 13 ))} )
test_that( "get_max_prob_future_map function ",  {expect_that(getBestPosition("AAAAAAAAAAGDS(0.206)DS(0.794)WDADTFSMEDPVRK")  ,
                                                              equals( 15))} )




##--------

context("get_pos_string")

a <- list(140,166,179	)
b <- list( 5,9)

test_that( "get_pos_string function ",  {expect_that(getPosString(a, b),
                                                     equals( "(144;148)|(170;174)|(183;187)" ))} )

a <- list(140	)
b <- list( 5,9)

test_that( "get_pos_string function ",  {expect_that(getPosString(a, b),
                                                     equals( "144;148" ))} )


##--------
context("get_pos_string")


aa_seq_fasta <- "MAALVLEDGSVLQGRPFGAAVSTAGEVVFQTGMVGYPEALTDPSYKAQILVLTYPLIGNYGIPSDEEDEFGLSKWFESSEIHVAGLVVGECCPTPSHWSANCTLHEWLQQRGIPGLQGVDTRELTKKLREQGSLLGKLVQKGTEPSALPFVDPNARPLAPEVSIKTPRVFNAGGAPRICALDCGLKYNQIRCLCQLGAEVTVVPWDHELDSQKYDGLFLSNGPGDPASYPGVVSTLSRVLSEPNPRPVFGICLGHQLLALAIGAKTYKMRYGNRGHNQPCLLVGTGRCFLTSQNHGFAVDADSLPAGWAPLFTNANDCSNEGIVHDSLPFFSVQFHPEHRAGPSDMELLFDVFLETVREAAAGNIGGQTVRERLAQRLCPPELPIPGSGLPPPRKVLILGSGGLSIGQAGEFDYSGSQAIKALKEENIQTLLINPNIATVQTSQGLADKVYFLPITLHYVTQVIRNERPDGVLLTFGGQTALNCGVELTKAGVLARYGVRVLGTPVETIELTEDRRAFAARMAEIGEHVAPSEAANSLEQAQAAAERLGYPVLVRAAFALGGLGSGFASTKEELSALVAPAFAHTSQVLIDKSLKGWKEIEYEVVRDAYGNCVTVCNMENLDPLGIHTGESIVVAPSQTLNDREYQLLRRTAIKVTQHLGIVGECNVQYALNPESEQYYIIEVNARLSRSSALASKATGYPLAYVAAKLALGIPLPELRNSVTGGTAAFEPSLDYCVVKIPRWDLSKFLRVSTKIGSCMKSVGEVMGIGRSFEEAFQKALRMVDENCVGFDHTVKPVSDMELETPTDKRIFVVAAALWAGYSVERLYELTRIDCWFLHRMKRIVTHAQLLEQHRGQALPQDLLHQAKCLGFSDKQIALAVLSTELAVRKLRQELGICPAVKQIDTVAAEWPAQTNYLYLTYWGNTHDLDFRAPHVLVLGSGVYRIGSSVEFDWCAVGCIQQLRKMGYKTIMVNYNPETVSTDYDMCDRLYFDEISFEVVMDIYELENPEGVILSMGGQLPNNMAMALHRQQCRVLGTSPEAIDSAENRFKFSRLLDTIGISQPQWRELSDLESARQFCHTVGYPCVVRPSYVLSGAAMNVAYTDGDLERFLSSAAAVSKEHPVVISKFIQEAKEIDVDAVACDGIVSAIAISEHVENAGVHSGDATLVTPPQDITPKTLERIKAIVHAVGQELQVTGPFNLQLIAKDDQLKVIECNVRVSRSFPFVSKTLGVDLVALATRIIMGEKVEPVGLMTGSGVVGVKVPQFSFSRLAGADVVLGVEMTSTGEVAGFGESRCEAYLKAMLSTGFKIPEKNILLTIGSYKNKSELLPTVRLLESLGYSLYASLGTADFYTEHGVKVTAVDWHFEEAVDGECPPQRSILDQLAENHFELVINLSMRGAGGRRLSSFVTKGYRTRRLAADFSVPLIIDIKCTKLFVEALGQIGPAPPLKVHVDCMTSQKLVRLPGLIDVHVHLREPGGTHKEDFASGTAAALAGGVTMVCAMPNTRPPIIDAPALALAQKLAEAGARCDFTLFLGASSENAGTLGAVAGSAAGLKLYLNETFSELRLDSVAQWMEHFETWPAHLPIVAHAERQSVAAVLMVAQLTQRPVHICHVARKEEILLIKTAKAQGLPVTCEVAPHHLFLNREDLERLGPGKGEVRPELGSREDMEALWENMAVIDCFASDHAPHTLEEKCGPKPPPGFPGLETMLPLLLTAVSEGRLSLDDLLQRLHHNPRRIFHLPLQEDTYVEVDLEHEWTVPSHMPFSKARWTPFEGQKVKGTVRRVVLRGEVAYIDGQVLVPPGYGQDVRKWPQGVVPQPPPSTPATTEITTTPERPRRVIPGLPDGRFHLPPRIHRASDPGLPAEEPKEKPPRKVVEPELMGTPDGPCYPAPPVPRQASPQNLGSSGLLHPQMSPLLHSLVGQHILSVKQFTKDQMSHLFNVAHTLRMMVQKERSLDILKGKVMASMFYEVSTRTSSSFAAAMARLGGAVLSFSEATSSVQKGESLADSVQTMSCYADVIVLRHPQPGAVELAAKHCRRPVINAGDGVGEHPTQALLDIFTIREELGTVNGMTITMVGDLKHGRTVHSLACLLTQYRVSLRYVAPPSLRMPPSVRDFVASRGTKQEEFESIEEALPDTDVLYMTRIQKERFGSVQEYEACFGQFILTPHIMTRAKKKMVVMHPMPRVNEISVEVDSDPRAAYFRQAENGMYIRMALLATVLGRF"

list_of_15_mers <- c( "MAALVLEDGSVLQGR",
                      "_MAALVLEDGSVLQG",
                      "__MAALVLEDGSVLQ",
                      "___MAALVLEDGSVL",
                      "____MAALVLEDGSV",
                      "_____MAALVLEDGS",
                      "______MAALVLEDG",
                      "LATVLGRF_______",
                      "MYIRMALLATVLGRF",
                      "GMYIRMALLATVLGR" )

positions <- c( 8:2, 2225, 2218, 2217  )

for ( i in seq_along(positions) ) {

  position <- positions[i]
  my_15_mer <- list_of_15_mers[i]

  test_that( "get_pos_string function ",  {expect_that(getXMerString(aa_seq_fasta, "B2RQC6", position ),
                                                       equals( my_15_mer ))} )
}

#purrr::walk(  positions, list_of_15_mers,  test_function  )


##--------
context("formatPhosphositePosition")

test_that ( "Test formatPhosphositePosition 1", {expect_that(formatPhosphositePosition( "MFG1", "(1986)|(1998)|(2010)	", "S"), equals("MFG1:S1986")) })
test_that (  "Test formatPhosphositePosition 2", {expect_that(formatPhosphositePosition( "MFG1", "1082;1087", "S;Y"), equals("MFG1:S1082,Y1087")) })

##--------
