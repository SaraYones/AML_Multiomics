rules = (read.csv2('/Users/sarayones/Desktop/NetworksAllGenes-2019-01-07.txt', sep='\t', header = FALSE, col.names = c('FEATURES', 'DECISION', 'ACC_RHS', 'SUPP_RHS'), stringsAsFactors=FALSE))

rules$ACC_RHS = as.numeric(rules$ACC_RHS)
rules$SUPP_RHS = as.numeric(rules$SUPP_RHS)
#in this moment adding p-value column is necessary. You can add random values.
rules$PVAL = 0.05
visunet_out2 = visunet(rules, 'L')


out = rosetta(autcon)
rules = out$main
visunet_out = visunet(rules)