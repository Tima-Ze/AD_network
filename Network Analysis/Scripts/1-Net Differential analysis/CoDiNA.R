library(dplyr)
#set a directory to save the output
output_dir <- 'Results/codina/'

if(!dir.exists(file.path(output_dir))) {
  dir.create(file.path(output_dir),recursive = TRUE)}

################ Data preparation ###################

#read the consensus networks
for (i in list.files(path='Data/Networks/', pattern='*txt', full.names = T)){
  if (grepl('AD', i)){
    assign('ad', 
           data.table::fread(i) %>% 
             dplyr::rename(wTO=CN))
  }
  else{
    assign('ctrl', 
           data.table::fread(i) %>% 
             dplyr::rename(wTO=CN))
  }
}

#Since we want CoDiNA to keep all links into comparison analysis,
#we need to add the missing nodes to the networks with assigning a weight of 0
#This is because CoDiNA will remove the nodes that are not present in both networks

node_ctrl <- unique(c(ctrl$Node.1, ctrl$Node.2))
node_ad <- unique(c(ad$Node.1, ad$Node.2))

add_ctrl <- setdiff(node_ad, node_ctrl)
add_ad <- setdiff(node_ctrl,node_ad)

df_ad <- data.frame(Node.1=add_ad, Node.2=add_ad, wTO=0)
df_ctrl <- data.frame(Node.1=add_ctrl, Node.2=add_ctrl, wTO=0)

#Make sure there are not duplicated links
ctrl <- rbind(ctrl, df_ctrl) %>% 
  distinct()
ad <- rbind(ad,df_ad)%>% 
  distinct()


################ CoDiNA ###################
library(CoDiNA)
DiffNet = MakeDiffNet(Data = list(ad, ctrl),
                      Code = c("AD", 'Ctrl'),
                      stretch=F)

write.table(DiffNet, paste0(output_dir,'DiffNet.txt'), sep = '\t', quote = F, row.names = F)

################ Nodes ###################
#Using median:
int_C = quantile(DiffNet$Score_internal, 0.5)
ext_C = quantile(DiffNet$Score_Phi, 0.5)

Nodes_Groups = ClusterNodes(DiffNet = DiffNet, 
                            cutoff.external = ext_C, 
                            cutoff.internal = int_C)
table(Nodes_Groups$Phi_tilde)

write.table(Nodes_Groups, paste0(output_dir,"DiffNode_median.txt"), sep="\t", quote = F)
