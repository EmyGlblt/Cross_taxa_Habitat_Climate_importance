library(sjPlot)

colnamesD = c('Habitat type',	'CLC categories',	
             'Proportion \n within study sites',	'Proportion \n within Finland')

col1 = c('Urban',
         'Agricultural areas',
        ' Semi natural herbaceous',
         'Forest',
         'Wetlands',
         'Water bodies')

col2 = c('111 to 142',
         '210, 220, 241, 242',
         '231, 243, 244, 321, 322, \n 324, 331, 332, 333',
         '311, 312, 313',
         '411 to 423',
         '511 to 523')

col3 = c(18.09,
         10.12,
         12.17,
         41.64,
         16.44,
         1.54)

col4 = c(2.75,
         5.99,
         14.26,
         48.39,
         6.68,
         21.9)


Df = as.data.frame(cbind(col1, col2, col3, col4))
colnames(Df) = colnamesD


sjtab(Df) 


library("rempsyc")

nice_table(Df)
