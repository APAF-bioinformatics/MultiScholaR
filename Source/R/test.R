
# Test voom-limma
```{r}
p_load(edgeR)
p_load(RUVSeq)

#filtered_norm <- calcNormFactors(filtered_mat, method = "TMM")


ruv_mat <- design_matrix %>% 
  group_by( group) %>%
  mutate( replicate_id = row_number()) %>%
  ungroup() %>% 
  pivot_wider( id_cols = group,
               names_from = replicate_id,
               values_from = Sample_ID) %>%
  column_to_rownames("group") %>%
  as.matrix()

seqRUVs <- RUVs(log2(as.matrix(filtered_mat)), control_genes_index, k=5, ruv_mat, isLog=TRUE)


filtered_voom <- voom( filtered_mat, cbind( design_m, seqRUVs$W ), plot=TRUE, normalize.method = "cyclicloess")


cbind( design_m, seqRUVs$W )

vfit <- lmFit(filtered_voom, cbind( design_m, seqRUVs$W ))

contr.matrix.ruv <- makeContrasts(contrasts = contrast_strings,
                              levels = colnames(cbind( design_m, seqRUVs$W )))


vfit <- contrasts.fit(vfit, contrasts=contr.matrix.ruv)
efit <- eBayes(vfit)

efit

plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))


de_tbl_test <- topTreat(efit, coef = contrast_strings[[1]], n = Inf) %>%
  mutate( q.mod= qvalue(P.Value)$q) %>%
  mutate( fdr.mod  = p.adjust(P.Value, method="BH")) %>%
  dplyr::rename( p.mod = P.Value)

voom_tbl <- de_tbl_test %>%
  rownames_to_column("uniprot_acc") %>%
  dplyr::select( uniprot_acc, logFC, q.mod)


orig_tbl <- de_tbl %>%
  rownames_to_column("uniprot_acc") %>%
  dplyr::select( uniprot_acc, logFC, q.mod)



compare_orig_vroom <- voom_tbl %>%
  left_join( orig_tbl, by = c("uniprot_acc"), suffix=c(".orig", ".voom"))

compare_orig_vroom %>%
  ggplot(aes( logFC.orig, logFC.voom)) +
  geom_point()

compare_orig_vroom %>%
  ggplot(aes( q.mod.orig, q.mod.voom)) +
  geom_point()

```