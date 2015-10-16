package GoFactory;


use strict;
use warnings;
use vars qw($AUTOLOAD);
use DBI;
use GO::AppHandle;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(createHandler getPaths getlistterm getUp getAllup getDowntodepth getParents);
our $VERSION = '1.00';



my $apph=createHandler();

#--------------------------------
#fct createHandler()
#return un apphandler
#--------------------------------
sub createHandler {

#  my $dbname = "go";
#  my $user="flegeai";
#  my $password = "toto";
#  my $apph = GO::AppHandle->connect(-dbname=>$dbname,-dbuser=>$user, -dbauth=>$password);

  # Connection at GO
#  my $dbname = "go";
#  my $mysqlhost = "pythagore";
#  my $dbuser ='sgo';
#  my $dbpwd = 'selGo2005';
#  my $port = 3306;

  my $dbname = "go0207";
  my $mysqlhost = "pythagore.versailles.inra.fr";
  my $dbuser ="sgo0207";
  my $dbpwd = "sgo0207";
  my $port = 3306;


  my $apph = GO::AppHandle->connect(-dbname=>$dbname, -dbhost=>$mysqlhost, -dbuser=>$dbuser, -dbauth=>$dbpwd, -port=>$port);


  $apph;

}


#--------------------------------
#fct getPaths(term)
#return refarray : paths for term
#--------------------------------
sub getPaths{
  my $termref=shift @_;

  my $paths;

  if($termref->acc ne "GO:0003673"){

    my $graph=$apph->get_graph_by_terms([$termref], my $depth);
    $paths = $graph->paths_to_top($termref->acc);  
  }

  $paths;

}

#--------------------------------------------------------
#genere des liste de termacc et termname 
#a partir de tous les pathUp d'ter
#return %hgoid=(acc=> @tabtermaccs, name=>@tabtermnames)
#--------------------------------------------------------
sub getlistterm{

  my $termref=shift @_;

  my @tablistterms;
  my $paths=getPaths($termref);

  foreach my $path (@{$paths}){
    my $listterm=$path->term_list;
    unshift @{$listterm},($termref);
    my $lastterm=pop @{$listterm};
    if($lastterm->acc ne 'all'){print STDERR "\nATTENTION la racine est ",$lastterm->acc,"\n";}
    push @tablistterms,($listterm);
  }
  my @tabtermaccs;
  my @tabtermnames;
  my @tabterms;
  my %hgoid=(acc=>\@tabtermaccs, 
	    name=>\@tabtermnames,
	    term=>\@tabterms);

  foreach my $listterms (@tablistterms){
    my @revtabterm=reverse @{$listterms};
    foreach my $term (@revtabterm){ 
      my $termacc=$term->acc;
      my $termname=$term->name;
      if(grep(/$termacc/,@tabtermaccs) == 0){
	push @tabtermaccs,($termacc);
	push @tabtermnames,($termname);
	push @tabterms,($term);
      }
    } #term suivant
  } #listterm suivant
  %hgoid;
}


#--------------------------------------------------------
#imprime tous les chemin compactes du graph du termref vers le haut 
#--------------------------------------------------------
sub getUp{
  my $termref=shift @_;
  my $evidence=shift @_;
  local *FH =shift @_;

  print FH "GO_ID=",$termref->acc,"\tEV_CODE=",$evidence,"\tGO_NAME=",$termref->name,"\tGO_TYPE=",$termref->type,"\n";
  if($termref->acc ne  'all'){

    my $depthref=0;
    my $graph=$apph->get_graph_by_terms([$termref], my $depth);
    #print "TERM_COUNT=",$graph->node_count,"\t","PARENT=",$graph->n_parents($termref->acc),"\n";


    my $it = $graph->create_iterator;  #returns a GO::Model::GraphIterator object
    $it->compact(1);
    while (my $ni = $it->next_node_instance) {
      
      $depth = $ni->depth;
      my $term = $ni->term;

      if($depthref==0 || $depth<=$depthref){
	printf FH
	  "%s Term = %s (%s)   // depth=%d\n",
	  "--" x $depth,
	  $term->name,
	  $term->public_acc,
	  $depth;
	if($term->acc eq $termref->acc){
	  $depthref=$depth;
	  next;
	}
      }      
      
    } # of while next_node_instance
  } # of if ne "GO:0003673"

}

#--------------------------------------------------------
#imprime tous les chemin du graph du termref vers le haut 
#--------------------------------------------------------
sub getAllup{
  my $termref=shift @_;
  

  print "GO_ID=",$termref->acc,"\tGO_NAME=",$termref->name,"\tGO_TYPE=",$termref->type,"\n";

  if($termref->acc ne "GO:0003673"){
    my $depthref=-1;
    my $graph=$apph->get_graph_by_terms([$termref],my $depth);
    
    my $it = $graph->create_iterator;
    # returns a GO::Model::GraphIterator object
    
    while (my $ni = $it->next_node_instance) {
      $depth = $ni->depth;
      my $term = $ni->term;
      #
      if($depthref==-1 || $depth<=$depthref){
	#if($depth<=3 || $depth==$depthref){
	printf
	  "%s Term = %s (%s)   // depth=%d\n",
	  "--" x $depth,
	  $term->name,
	  $term->public_acc,
	  $depth;
	if($term->acc eq $termref->acc){
	  $depthref=$depth;
	  #print $depthref,"\n";
	  next;
	}
	#}
      }      
      
    }
  }

}

#--------------------------------------------------------
#imprime tous les chemin du graph du termref jusqu'a depth 
#--------------------------------------------------------
sub getDowntodepth{
  my $termref=shift @_;
  local *FH =shift @_;
  my $depthref=shift @_;

  #print "GO_ID=",$termref->acc,"\tGO_NAME=",$termref->name,"\tGO_TYPE=",$termref->type,"\n";

  if($termref->acc ne 'all'){
    my $graph=$apph->get_graph_by_terms([$termref],my $depth);
    
    my $it = $graph->create_iterator;
    # returns a GO::Model::GraphIterator object
    
    while (my $ni = $it->next_node_instance) {
      $depth = $ni->depth;
      my $term = $ni->term;
      #
      if($depth<=$depthref){
	printf
	  "%s Term = %s (%s)   // depth=%d\n",
	  "--" x $depth,
	  $term->name,
	  $term->public_acc,
	  $depth;
      } else{
	next;
      }
      
    }
  }

}
#-------------------------------------------------
#imprime 1 seul chemin via get_parent_terms(termi)
#--------------------------------------------------
sub getParents {
  
  my $term=shift @_;
  my @tabparent;
  
  my $term_lref= $apph->get_parent_terms($term);
  push  @tabparent, ($term_lref);
  
  while(my $term_lref->[0] != "") { 
    if($term_lref->[0] != ""){
      $term_lref = $apph->get_parent_terms($term_lref->[0]);
      if($term_lref->[0] != "") {
	push @tabparent,($term_lref);
      }
    }
  }
  
  my $depth=0;

  for(my $p=scalar(@tabparent)-2; $p>=0;$p--) {
    my $terms=$tabparent[$p];  
    if(scalar(@{$terms})){ 
      foreach my $term (@{$terms}){
	print "----" x $depth, $term->acc,"\t",$term->name,"\n";
      }
    }
    $depth++;
  }  
}

