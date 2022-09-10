


"uniprot_acc"	"reactome_id"	"pathway_url"	"pathway_name"	"evidence_type"	"species"




reactome_file<-file.path(args$tmp_dir,args$reactome_file)
if(!file.exists(reactome_file))
{
  logwarn("Download Reactome UniProt to pathways file.")
  status <- download.file(url="https://reactome.org/download/current/UniProt2Reactome.txt", destfile=reactome_file)
  loginfo(status)
}
loginfo("Reading Reactome UniProt to pathways file.")
captured_output<-capture.output(
  reactome_map <- vroom::vroom( reactome_file ,
                                col_names = c("uniprot_acc", "reactome_id", "url", "reactome_term", "evidence", "organism") )
  , type = "message"
)
logdebug(captured_output)

