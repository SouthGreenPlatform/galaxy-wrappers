// ------------------------------------------------------------------
// Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
//
// This program is open source; you can redistribute it and/or modify
// it under the terms of the Artistic License (see LICENSE file).
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
//
// You should have received a copy of Artistic License along with
// this program; if not, please see http://www.opensource.org
//
// $Id: Param.cc,v 1.54 2009-01-12 14:07:37 sallet Exp $
// ------------------------------------------------------------------
// File:     Param.cc
// Contents: class Parameters
// ------------------------------------------------------------------

#include "Param.h"


// -----------------------------------------------
//  Gestion des arguments (positive integer)
// -----------------------------------------------
int TestIArg(char *arg)
{
  int tmp;
  return (sscanf(arg, "%d", &tmp) == 1) && (tmp >= 0);
}
// -----------------------------------------------
//  Gestion des arguments (arbitrary integer)
// -----------------------------------------------
int TestIAnyArg(char *arg)
{
  int tmp;
  return (sscanf(arg, "%d", &tmp) == 1);
}

// -----------------------------------------------
//  Gestion des arguments (double)
// -----------------------------------------------
int TestDArg(char *arg)
{
    double tmp;
    return  (sscanf(arg, "%lf", &tmp) == 1) && (tmp >= 0);
}

// ------------------------
//  initParam.
// ------------------------
void Parameters :: initParam (int argc, char * argv[])
{
  char *key, *val = NULL;
  key = new char[FILENAME_MAX+1];
  val = new char[FILENAME_MAX+1];
  char *key2, *val2 = NULL;
  key2 = new char[FILENAME_MAX+1];
  val2 = new char[FILENAME_MAX+1];

  fprintf(stderr,"EuGene rel. %s ",VERSION);
  fflush(stderr);

  // store in the map the eugene directory
  strcpy(key,"eugene_dir");
  if (EUGENE_DIR == NULL) 
    m["eugene_dir"] = DEFAULT_EUGENE_DIR;
  else 
    m["eugene_dir"] = EUGENE_DIR;
  strcpy(val, getC("eugene_dir"));
  strcat(val,"/");
  m[key] = val;

  UpdateParametersFileName(argc, argv);
  ReadPar((std::string) getC("parameters_file"));
  ReadArg(argc, argv);

  if (m.count("Param.debug")) 
    for (iter = m.begin(); iter!=m.end(); ++iter)
      fprintf(stderr,"%s = %s\n",iter->first, iter->second);

  iter = m.begin();  // Cf. : getUseSensor
  
  fprintf(stderr,"- %s -\n",getC("EuGene.organism"));
  fprintf(stderr,"Parameters file %s loaded.\n", getC("parameters_file"));
  if (getI("EuGene.sloppy")) fprintf(stderr,"WARNING: Sloppy mode selected.\n");
  strcpy(key2,"web_dir");
  if (!strcmp(getC("Output.webdir"),"LOCAL")) {
    strcpy(val2, getC("eugene_dir"));
    strcat(val2,WEB_DIR);
  } else 
    strcpy(val2, getC("Output.webdir"));
  m[key2] = val2;
}

// ---------------------------------------------------
// check if a parameters file name is set in arguments
// ---------------------------------------------------
void Parameters :: UpdateParametersFileName(int argc, char * argv[])
{
  int carg;
  char *key, *val;
  key = new char[FILENAME_MAX+1];
  val = new char[FILENAME_MAX+1];
  std::string name = "";

  while ((carg = getopt(argc, argv, POSSIBLE_ARGUMENTS)) != EOF) {
    switch (carg) {
    case 'A': 
      if (optarg) name = optarg;
      else {ShowUsage();  exit(1);}
    }
  }
  
  if (name == "") 
    name = ((std::string) getC("eugene_dir"))+DEFAULT_PARA_FILE;

  strcpy(key,"parameters_file");
  strcpy(val, name.c_str());
  m[key] = val;
}

// ------------------------
//  Read arguments line.
// ------------------------
void Parameters :: ReadArg(int argc, char * argv[])
{
  int carg, errflag = 0;
  char *key, *val = NULL;
  char* indexPos = NULL;

  m["fstname"] = "\000"; // no default input

  optind = 1;            // reinit this getopt static variable
 
  while ((carg = getopt(argc, argv, POSSIBLE_ARGUMENTS)) != EOF) {
    switch (carg) {
      
    case 'D':           /* Definition of any parameter */

      if ( (indexPos = index(optarg,'=')) && (indexPos > optarg) && (indexPos < optarg+strlen(optarg)-1)) {
	key = new char[FILENAME_MAX+1];
	val = new char[FILENAME_MAX+1];
	*indexPos = 0;
	sscanf(optarg, "%s", key);
	sscanf(indexPos+1, "%s", val);
	m[key] = val;
      }
      else {
	fprintf(stderr,"\nInvalid command line parameter definition: %s\n",optarg);
	exit(2);
      }
      break;

	case 'F': /* Remove fragment proteins from output */
	m["Output.RemoveFrags"] = "1";
	break;

    case 'n':           /* -n normalize across frames */
      if (!TestIArg(optarg))
	errflag++;
      else m["Output.normopt"] = optarg;
      break;
      
    case 'h':           /* help */
      errflag++;
      break;
      
    case 's':           /* Single gene mode. Start/End on IG only */
      m["EuGene.ExonPrior"]       = "0.0";
      m["EuGene.IntronPrior"]     = "0.0";
      m["EuGene.FivePrimePrior"]  = "0.0";
      m["EuGene.ThreePrimePrior"] = "0.0";
      break;

    case 'r':    /* RepeatMasker input */
      m["Sensor.Repeat.use"] = "1";
      break;
      
    case 'a':    /* AltEst activation */
      m["AltEst.use"] = "1";
      break;
      
    case 'c':           /* -c couverture */
      if (! TestIArg(optarg))
        errflag++;
      else m["Output.golap"] = optarg;
      break;
      
    case 'l':           /* -l imglen */
      if (! TestIArg(optarg))
        errflag++;
      else m["Output.glen"] = optarg;
      break;
      
    case 'u':           /* -u From */
      if (! TestIArg(optarg))
	errflag++;
      else m["Output.gfrom"] = optarg;
      break;
      
    case 'v':           /* -v To */
      if (! TestIArg(optarg))
	errflag++;
      else m["Output.gto"] = optarg;
      break;
      
    case 'U':
      if (! TestIArg(optarg)) 
	errflag++;
      else {
        int tmp;
        key = new char[FILENAME_MAX+1];
        sscanf(optarg, "%d", &tmp);
	sprintf(key,"%d",Max(0,tmp-1));
	m["EuGene.from"] = key;
	}
      break;

    case 'V':
      if (! TestIArg(optarg)) 
	errflag++;
      else {
        int tmp;
        key = new char[FILENAME_MAX+1];
        sscanf(optarg, "%d", &tmp);
	sprintf(key,"%d",Max(0,tmp-1));
        m["EuGene.to"] = key;
	}
      break;

    case 'x':           /* -x resx */
      if (! TestIArg(optarg))
        errflag++;
      else m["Output.resx"] = optarg;
      break;
      
    case 'y':           /* -y resy */
      if (! TestIArg(optarg))
        errflag++;
      else m["Output.resy"] = optarg;
      break;
      
    case 'g':          /* -g activate graphical output */
      m["Output.graph"] = "1";
      break;

    case 'O':           /* -O output */
      m["Output.Prefix"] = optarg;
      break;
      
    case 'p':           /* -p print opt: s -> short, l -> long,   d -> detailed
			                 g -> GFF,   a -> araset, o -> stdout(l) 
                                         h -> html */

      m["Output.format"] = optarg;
      for (int i=0; i<strlen(optarg); i++) {
	if (getI("Output.graph") == 0  &&  (optarg[i] == 'h')) {
	  m["Output.graph"] = "1"; // HTML output means graphical output 
	}
	if ((optarg[i] != 's') && (optarg[i] != 'l') && (optarg[i] != 'g') &&
	    (optarg[i] != 'd') && (optarg[i] != 'h') && (optarg[i] != 'a') &&
	    (optarg[i] != 'o'))
	  errflag++;
      }
      break;
      
    case 'm':           /* -m IMM matrix */
      m["Sensor.MarkovIMM.use"] = "1" ;
      m["MarkovIMM.matname"] = optarg;
      break;

    case 'M':           /* -M proteic markovian matrix */
      m["Sensor.MarkovProt.use"] = "1";
      m["MarkovProt.matname"] = optarg;
      break;

    case 'o':           /* -o offset */
      if (! TestIAnyArg(optarg))
	errflag++;
      else m["Output.offset"] = optarg;
      break;
      
    case 'w':           /* -w window */
      if (! TestIArg(optarg))
	errflag++;
      else m["Output.window"] = optarg;
      break;

    case 'b':           /* -b  use blastx result */
      m["Sensor.BlastX.use"] = "1";
      if (optarg) {
	m["BlastX.levels"] = optarg;
	if(strlen(optarg) > 10)
	  errflag++;
      }
      break;

    case 'E':
      m["Est.PostProcess"] = "2";
      break;

    case 'B':
      m["BlastX.PostProcess"] = "2";
      break;

    case 'G':           /* -G use sensor GFF */
      m["Sensor.GFF.use"] = "1";
      break;

    case 'd':           /* -d  use cDNA blastn results      */
      m["Sensor.Est.use"] = "1";
      break; 

    case 'R':           /* -R use RAFL-like EST */
      m["Sensor.Riken.use"] = "1";
      break;

    case 'f':           /* -f frameshift loglike */
      if (! TestDArg(optarg))
	errflag++;
      else m["EuGene.frameshift"] = optarg;

      break;
      
    case 't':           /* -t use tblastx results */
      m["Sensor.Homology.use"] = "1";
      if (optarg) {
	m["Homology.protmatname"] = optarg;
      }
      break;

    case 'Z':           /* parameters optimization */
      m["ParaOptimization.Use"] = "1";
      if (optarg) m["ParaOptimization.TrueCoordFile"] = optarg;
      break;

    case '?':           /* bad option */
      errflag++;
    }
  }

  // may remain arguments -> fasta filenames
  if ((argc - optind) == 0)
    errflag++;
  
  // check usage
  if (errflag) { ShowUsage();  exit(1); }
  
  // Set number of sequences given in arguments in the map
  setD("NbSequence", argc - optind); 
}

// ------------------------
// Show usage
// ------------------------
void Parameters :: ShowUsage (void)
{
  fprintf(stderr, "\nUsage: EuGene [options] FASTA-files\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  General:\n");
  fprintf(stderr, "    -A paramfile            Selects another parameter file\n");
  fprintf(stderr, "    -a                      Activates the alternative splicing prediction using EST\n");
  fprintf(stderr, "    -D<parameter>=<value>   Changes the parameter value\n"); 
  fprintf(stderr, "    -s                      Forbids partial gene prediction\n");
  fprintf(stderr, "    -U                      Predict only from this nucleotide position\n");
  fprintf(stderr, "    -V                      Predict only up to this nucleotide position\n");
  fprintf(stderr, "    -Z                      Select the parameters optimization mode\n");
  fprintf(stderr, "    -h                      Print this help\n\n");

  fprintf(stderr, "  Output control:\n");
  fprintf(stderr, "    -p a|d|g|h|l|s|o        Changes text output format (-pgh: html and gff)\n");
  fprintf(stderr, "    -g                      Activates graphical output\n");
  fprintf(stderr, "    -c overlap              Controls overlap between images\n");
  fprintf(stderr, "    -l length               Length of sequence per image\n");
  fprintf(stderr, "    -u position             Specifies the first base plotted\n");
  fprintf(stderr, "    -v position             Specifies the last base plotted\n");
  fprintf(stderr, "    -x resolution           Specifies the horizontal resolution\n");
  fprintf(stderr, "    -y resolution           Specifies the vertical resolution\n");
  fprintf(stderr, "    -w width                Changes the score smoothing window width\n");
  fprintf(stderr, "    -n 0|1|2                Changes the normalization of scores in graphs\n");
  fprintf(stderr, "    -o offset               Offset the base count in text/images\n");
  fprintf(stderr, "    -O dir                  Specifies an output directory\n\n");

  fprintf(stderr, "  Plugin control:\n");
  fprintf(stderr, "    -b levels               Activates BlastX plugin\n");
  fprintf(stderr, "    -B                      Activates postprocessing in BlastX plugin\n");
  fprintf(stderr, "    -d                      Activates the EST plugin\n");
  fprintf(stderr, "    -E                      Activates postprocessing in EST plugin\n");
  fprintf(stderr, "    -f                      Changes the frameshift penalty\n");
  fprintf(stderr, "    -G                      Activates the GFF plugin\n");
  fprintf(stderr, "    -m DNA Markov matrix    Changes the Markov matrix used (IMM)\n");
  fprintf(stderr, "    -M AA  Markov matrix    Activates the MarkovProt plugin\n");
  fprintf(stderr, "    -r                      Activates the Repeat plugin\n");
  fprintf(stderr, "    -R                      Activates the Riken (FL cDNA) plugin\n");
  fprintf(stderr, "    -t AA similarity matrix Activates the Homology plugin\n\n");
	
  exit(1);
}

// ------------------------
//  Read parameters file.
// ------------------------ 
void Parameters :: ReadPar(std::string  para_file_name)
{
  char line    [MAX_LINE];
  char *key, *val = NULL;
  int  n=0;

  fp = FileOpen(NULL,para_file_name.c_str(),"r");
  
  while(fgets (line, MAX_LINE, fp) != NULL) {
    n++;
    if (line[0] != '#') {
      key = new char[FILENAME_MAX+1];
      val = new char[FILENAME_MAX+1];
      if(sscanf(line, "%s %s", key, val) == 2) {
        // SMachine.cmd		"splicemachine.pl "
	if (val[0] == '"') {
	  char *from = index(line,'"');
	  char *to = rindex(from,'"');
	  if (from == to) {
	    fprintf(stderr, "\nIncorrect parameter file %s line %d\n",
		    para_file_name.c_str(),n);
	    exit(2);
	  }
	  *to = 0;
	  strcpy(val,from+1);
	}
	m[key] = val;
      }
      else {
	fprintf(stderr, "\nIncorrect parameter file %s line %d\n",
		para_file_name.c_str(),n);
	exit(2);
      }
    }
  }

  fclose(fp);
  // Check eugene version and the parameter file version are the same 
  if(strcmp(getC("EuGene.version"), VERSION)) {
    fprintf(stderr, "\nIncorrect parameter file version : %s\n", getC("EuGene.version"));
    fprintf(stderr,"Version %s required\n", VERSION);
    exit(2);
  }
}

// ------------------------
//  count.
// ------------------------
int Parameters :: count(char *key)
{
  return ( m.count(key) );
}

// ------------------------
//  test if a parameter is 
//  defined w/o issuing a warning
// ------------------------
bool Parameters::probeKey(char *key, int index){

  if (!index) return m.count(key);

  int len = strlen(key);
  char *altkey = new char[len+10];

  strcpy(altkey,key);
  key = altkey+len;
  sprintf(key,"[%d]",index);

  if (m.count(altkey)) {
    delete [] altkey; 
    return true;
  }
  return false;
}

// ------------------------
//  getChar param.
// ------------------------
char* Parameters :: getC(const char *key, int index, int sloppy)
{
  // case where there is a unique value for 'key' 
  if (!index && m.count(key))
    return (char*)m[key];

  int len = strlen(key);
  char *altkey = new char[len+10];
  char * charptr;

  strcpy(altkey,key);
  sprintf(altkey+len,"[%d]",index);
  // get the value for key[index]
  if (m.count(altkey)) {
    charptr = (char*)m[altkey];
    delete [] altkey;

    return charptr;
  }
  else 
  {
    if (sloppy) { return 0; }
    else {
      fprintf(stderr,"\nError: Undefined key %s\n",altkey);
      exit(2);
    }
  }
}

// ------------------------
//  getDouble param.
// ------------------------
double Parameters :: getD(const char *key, int index, int sloppy)
{
  char *res = getC(key,index,sloppy);
  double dou = 0.0;
  if (res) 
  {
    if(!strcmp(res, "NINFINITY"))
      dou = NINFINITY;
    else
      dou = atof(res);
  }
  return dou;
}

// ------------------------
//  getInt param.
// ------------------------
int Parameters :: getI(const char *key, int index, int sloppy)
{
  char *res = getC(key,index,sloppy);
  
  if(res) {
    if(!strcmp(res, "TRUE"))
      return 1;
    else if(!strcmp(res, "FALSE"))
      return 0;
    else
      return atoi(res);
  }
  return 0;
}
  
// ------------------------
//  Get Use.Sensor.
//  Put in key and val respectively the name and the priority of the next used sensor
// ------------------------
int Parameters :: getUseSensor(char **key, int *val)
{
  int l;
  char s[FILENAME_MAX+1];

  while(iter != m.end()) {
    while(iter != m.end() &&
	  (iter->first[0] != 'S' || iter->first[1] != 'e' ||
	   iter->first[2] != 'n' || iter->first[3] != 's' ||
	   iter->first[4] != 'o' || iter->first[5] != 'r'))
      ++iter;
    if(iter != m.end()) {
      l = strlen(iter->first) - 1;
      if(iter->first[l-2] == 'u' && iter->first[l-1] == 's' && iter->first[l] == 'e')  
	++iter;
      else {
	strcpy(s,iter->first); strcat(s,".use"); 
	if (getI(s,0,getI("EuGene.sloppy"))) {
	  *key = (char*)iter->first;
	  *val = atoi(iter->second);
	  ++iter;
	  return 1;
	}
	++iter;
      }
    }
  }
  return 0;
}

// ------------------------
//  set param.
// ------------------------
void Parameters :: set (const char *key, const char *value)
{
  m[key] = value;
}

void Parameters :: setD (const char *key, double n)
{
  char *buffer = new char[FILENAME_MAX+1];
  sprintf (buffer, "%10f",n);
  delete [] m[key];
  m[key] = buffer;
}

void Parameters :: setI (const char *key, int n)
{
  char *buffer = new char[FILENAME_MAX+1];
  sprintf (buffer, "%d",n);
  delete [] m[key];
  m[key] = buffer;
}

// ------------------------
//  Default destructor.
// ------------------------
Parameters :: ~Parameters ()
{
// Need to be commented, if not segmentation fault
//   for (iter = m.begin(); iter!=m.end(); ++iter) { 
//     delete [] iter->second;
//     delete [] iter->first;
//  }
  m.clear();
}

// ----------------------------------------------------------
// WriteParam: write a new parameter file, named <parameters_file>.<date>.OPTI,
//             modifying the value of parameters in para_name
//             with the new value in para_val
// BEWARE: the comments after a parameter value change are not kept
//         in case of no comment the result of strchr(line,'#') was: (null).
// Evaluation: the name of the written parameter file 
// -----------------------------------------------------------
std::string Parameters::WriteParam (std::vector<std::string> para_name, 
				    std::vector<double> para_val)
{
  FILE *fp, *fp_opti;
  char line[MAX_LINE], new_line[MAX_LINE];
  char key[MAX_LINE], val[MAX_LINE];
  char filename[FILENAME_MAX+1];
  bool find_para;
  char *d = new char[MAX_LINE];
  std::string s; int i;

  strcpy(filename, getC("parameters_file"));
  fp = FileOpen(NULL,filename,"r");

  // write the new file in the current repertory
  s = (std::string) filename; i = s.rfind("/");
  if (i != std::string::npos) s = s.substr(i+1,s.length());
  strcpy(filename, s.c_str());
  strcat(filename, ".");
  GetStrDate(d);
  strcat(filename, d);
  strcat(filename, ".OPTI");
  fp_opti = FileOpen(NULL,filename,"w");
  
  while(fgets (line, MAX_LINE, fp) != NULL) {
    find_para = false;
    if (line[0] != '#') 
      if (sscanf(line, "%s %s", key, val)==2) 
	for (i=0; i<para_name.size(); i++) 
	  if ( para_name[i] == (std::string) key ) {
	    find_para = true;
	    sprintf(new_line,"%s\t%f\n", key, para_val[i]);
	    i = para_name.size();
	  }
    if (!find_para)
      for (int i=0; i<MAX_LINE; i++) new_line[i] = line[i];
    fprintf(fp_opti,"%s",new_line);
  }

  fclose(fp);
  fclose(fp_opti);

  delete [] d;

  return (std::string) filename;
}



// -----------------------------------------------------------
// -----------------------------------------------------------
void Parameters::ResetIter(void)
{
  iter = m.begin();
}
