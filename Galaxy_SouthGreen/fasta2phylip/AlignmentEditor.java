import java.net.*;
import java.io.*;
import java.lang.reflect.*;
import java.util.*;

/**
* Executable program to edit / clean alignments
* @author	Dufayard Jean-Francois
* @version	1.0
* date : 01/2010
*/
public class AlignmentEditor {

// *******************************************************************************
//ATTRIBUTS
// *******************************************************************************

boolean proteic=false;

static final int CLUSTAL=0;

static int nbNo = 0;
static int nbUnsuf = 0;
static int nbNoGamma = 0;
static int nbNoBeta = 0;
static int nbValide = 0;

// *******************************************************************************

// *******************************************************************************
//EXECUTABLE & CONSTRUCTOR
// *******************************************************************************
	public static void main(String args[]) {
		AlignmentEditor alignmentEditor= new AlignmentEditor(args);
	}

	public AlignmentEditor(String args[]) {
		Alignment alignment=null;
		File input=null;
		String inputString=null;
		String output=null;
		String reference=null;
		String gamma=null;
		String beta=null;
		String stat=null;
		String suffix=null;
		String alternate=null;
		boolean strict=false;
		int operation=0;
		int nbBig=0;
		for (int i=0;i<args.length;i=i+2) {
		//System.out.println(args[i]);
			if (args[i].equalsIgnoreCase("-input")) {
				inputString=args[i+1];
				input= new File(args[i+1]);
			}
			if (args[i].equalsIgnoreCase("-output")) {
				output= args[i+1];
				//operation=1;
			}
			if (args[i].equalsIgnoreCase("-reference")) {
				reference= args[i+1];
			}
			if (args[i].equalsIgnoreCase("-proteic")) {
				proteic= true;
			}
			if (args[i].equalsIgnoreCase("-getBig")) {
				nbBig= (new Integer(args[i+1])).intValue();
			}
			if (args[i].equalsIgnoreCase("-getGapProportions")) {
				operation= 2;
				i--;
			}
			if (args[i].equalsIgnoreCase("-clustal2phylip")) {
				operation= 3;
				i--;
			}
			if (args[i].equalsIgnoreCase("-fasta2phylip")) {
				operation= 6;
				i--;
			}
			if (args[i].equalsIgnoreCase("-stats")) {
				operation= 5;
				i--;
			}
			if (args[i].equalsIgnoreCase("-strict")) {
				strict=true;
				i--;
			}
			if (args[i].equalsIgnoreCase("-custom")) {
				operation= 4;
				gamma= args[i+1];
				beta= args[i+2];
				stat= args[i+3];
				i++;
				i++;
			}
			if (args[i].equalsIgnoreCase("-getTaxaCode")) {
				operation= 7;
				gamma= args[i+1];
				beta= args[i+2];
				i++;
			}
			if (args[i].equalsIgnoreCase("-addsuffix")) {
				operation= 8;
				suffix= args[i+1];
				i++;
			}
			if (args[i].equalsIgnoreCase("-toUpperCase")) {
				operation= 9;
				i++;
			}
			if (args[i].equalsIgnoreCase("-removeN")) {
				operation= 10;
				i++;
			}
			if (args[i].equalsIgnoreCase("-ventileAlternateTranscripts")) {
				alternate= args[i+1];
				operation= 11;
				i++;
			}
		}
		//System.out.println(proteic);
		if (operation==11) {
			try {
				Hashtable parser= new Hashtable();
				Vector vector= new Vector();
				int numCluster=1;
				BufferedWriter writeUni = new BufferedWriter(new FileWriter(new File(output)));


				Alignment ali= new Alignment(input);
				for (int i=0;i<ali.sequences.size();i++) {
					String name= (String)(ali.names.elementAt(i));
					String sequence= (String)(ali.sequences.elementAt(i));
					if (name.indexOf(".")==-1) {
						writeUni.write( ">" + name + "\n");
						writeUni.flush();
						writeUni.write(((String)(ali.sequences.elementAt(i))) + "\n");
						writeUni.flush();
					} else {
						String identifier= name.substring(0,name.lastIndexOf("."));
						if (parser.containsKey(identifier)) {
							Vector local= (Vector)(parser.get(identifier));
							local.addElement(sequence);
							parser.put(identifier,local);
						} else {
							Vector local= new Vector();
							local.addElement(sequence);
							parser.put(identifier,local);
							vector.addElement(identifier);

						}
					}

				}
				for (int i=0;i<vector.size();i++) {
					String identifier= (String)(vector.elementAt(i));
					Vector local= (Vector)(parser.get(identifier));

					if (local.size()>1) {
						BufferedWriter write = new BufferedWriter(new FileWriter(new File(alternate + "/alternate" + numCluster + ".fasta")));
						String longer= "";
						for (int j=0;j<local.size();j++) {
							write.write( ">" + identifier + "." + (j+1) + "\n");

							write.flush();
							write.write(((String)(local.elementAt(j))) + "\n");
							write.flush();
							if (longer.length()<((String)(local.elementAt(j))).length()) {
								longer=((String)(local.elementAt(j)));
							}
						}
						numCluster++;
						write.close();
						writeUni.write( ">" + identifier + "." + local.size() + "\n");

						writeUni.flush();
						writeUni.write(longer + "\n");
						writeUni.flush();
					} else {
						writeUni.write( ">" + identifier + "." + local.size() + "\n");

						writeUni.flush();
						writeUni.write(((String)(local.elementAt(0))) + "\n");
						writeUni.flush();
					}



				}
				writeUni.close();
			} catch(Exception e) {
				e.printStackTrace();
			}
		} else if (operation==10) {
			try {

				Alignment ali= new Alignment(input);
				BufferedWriter write = new BufferedWriter(new FileWriter(new File(output)));
				for (int i=0;i<ali.sequences.size();i++) {
					int indexOf=((String)(ali.names.elementAt(i))).indexOf(" ");
					if (indexOf!=-1) {
						write.write( ">" + ((String)(ali.names.elementAt(i))).substring(0,indexOf) + "\n");
					} else {
						write.write( ">" + ((String)(ali.names.elementAt(i))) + "\n");
					}
					write.flush();
					write.write(((String)(ali.sequences.elementAt(i))).replaceAll("N","") + "\n");
					write.flush();
				}
				write.close();
			} catch(Exception e) {
				e.printStackTrace();
			}
		} else if (operation==9) {
			try {

				Alignment ali= new Alignment(input);
				BufferedWriter write = new BufferedWriter(new FileWriter(new File(output)));
				for (int i=0;i<ali.sequences.size();i++) {
					int indexOf=((String)(ali.names.elementAt(i))).indexOf(" ");
					if (indexOf!=-1) {
						write.write( ">" + ((String)(ali.names.elementAt(i))).substring(0,indexOf) + "\n");
					} else {
						write.write( ">" + ((String)(ali.names.elementAt(i))) + "\n");
					}
					write.flush();
					write.write(((String)(ali.sequences.elementAt(i))).toUpperCase() + "\n");
					write.flush();
				}
				write.close();
			} catch(Exception e) {
				e.printStackTrace();
			}
		} else if (operation==8) {
			try {

				Alignment ali= new Alignment(input);
				BufferedWriter write = new BufferedWriter(new FileWriter(new File(output)));
				for (int i=0;i<ali.sequences.size();i++) {
					int indexOf=((String)(ali.names.elementAt(i))).indexOf(" ");
					if (indexOf!=-1) {
						write.write( ">" + ((String)(ali.names.elementAt(i))).substring(0,indexOf) + "_" + suffix + "\n");
					} else {
						write.write( ">" + ((String)(ali.names.elementAt(i))) + "_" + suffix + "\n");
					}
					write.flush();
					write.write(((String)(ali.sequences.elementAt(i))) + "\n");
					write.flush();
				}
				write.close();
			} catch(Exception e) {
				e.printStackTrace();
			}
		} else if (operation==7) {
			try {
				BufferedReader buf1= new BufferedReader(new FileReader(new File(gamma)));
				BufferedReader buf2= new BufferedReader(new FileReader(new File(beta)));
				Vector correspondance= new Vector();
				String s= buf2.readLine();
				while (s!=null) {
					correspondance.addElement(s);
					s= buf2.readLine();
				}
				String ref= buf1.readLine();
				while (ref!=null) {
					int x=0;
					boolean founded=false;
					while (!founded && x<correspondance.size()) {
						String line= (String)(correspondance.elementAt(x));
						if (line.indexOf(ref)!=-1) {
							founded=true;
						}
						x++;
					}


					ref= buf1.readLine();
				}

			} catch(Exception e) {
				e.printStackTrace();
			}
		} else if (operation==4) {
			String ref1= null;
			String ref2= null;
			Hashtable tab1= new Hashtable();
			Hashtable tab2= new Hashtable();
			try {
				BufferedReader buf1= new BufferedReader(new FileReader(new File(gamma)));
				BufferedReader buf2= new BufferedReader(new FileReader(new File(beta)));
				BufferedWriter writeStat= new BufferedWriter(new FileWriter(new File(stat)));

				ref1 = buf1.readLine();
				if (ref1.endsWith("*")) {
					ref1=ref1.substring(0,ref1.length()-1);
				}
				while (ref1!=null) {
					if (!tab1.containsKey(ref1)) {
						tab1.put(ref1, new Integer(0));
					}
					ref1 = buf1.readLine();
					if (ref1!=null && ref1.endsWith("*")) {
						ref1=ref1.substring(0,ref1.length()-1);
					}
				}
				ref2 = buf2.readLine();
				if (ref2.endsWith("*")) {
					ref2=ref2.substring(0,ref2.length()-1);
				}
				while (ref2!=null) {
					if (!tab2.containsKey(ref2)) {
						tab2.put(ref2, new Integer(0));
					}
					ref2 = buf2.readLine();
					if (ref2!=null && ref2.endsWith("*")) {
						ref2=ref2.substring(0,ref2.length()-1);
					}
				}
				buf1.close();
				buf2.close();
				File inputFiles[] = (new File(inputString)).listFiles();
				String trace= "0 / " + inputFiles.length + " file computed...";
				System.out.print(trace);
				for (int i=0;i<inputFiles.length;i++) {
					if (inputFiles[i].getName().endsWith(".aln")) {
						Alignment ali= new Alignment(CLUSTAL,inputFiles[i]);
						if (i%100==0) {
							for (int k=0;k<trace.length();k++) {
								System.out.print("\b");
							}
							trace= i + " / " + inputFiles.length + " files computed...";
							System.out.print(trace);
						}

						boolean isOk = ali.filter(tab1,tab2,writeStat);

						if (isOk) {

							BufferedWriter write= new BufferedWriter(new FileWriter(new File(output + (inputFiles[i].getName()).substring(0,inputFiles[i].getName().lastIndexOf('.')) + ".fst")));
							write.write(ali.toString() + "\n");
							write.flush();
							write.close();
						}
					}


				}
				writeStat.write("#datasets gamma<3 and beta<3 : " + nbNo + "\n#datasets gamma>2 and beta=0 : " + nbNoBeta + "\n#datasets without gamma=0 and beta>2 : " + nbNoGamma + "\n#datasets selected : " + nbValide + "\n");
				writeStat.flush();
				writeStat.close();

			} catch(Exception e) {
				e.printStackTrace();
			}


		} else if (operation==6) {
			if (input.isDirectory()) {
				File inputFiles[] = (new File(inputString)).listFiles();
				for (int i=0;i<inputFiles.length;i++) {
					Alignment ali= new Alignment(inputFiles[i]);
					//System.out.println(ali);

					try {
						if (reference!=null) {
							BufferedReader read= new BufferedReader(new FileReader(new File(reference.replace("IN",(inputFiles[i].getName()).substring(0,inputFiles[i].getName().lastIndexOf('.'))))));
							String refString= read.readLine();
							ali.cleanAlignment(refString,inputFiles[i]);
						}
						BufferedWriter write= new BufferedWriter(new FileWriter(new File(output + (inputFiles[i].getName()).substring(0,inputFiles[i].getName().lastIndexOf('.')) + ".phy")));
						write.write(ali.phylipString() + "\n");
						write.flush();
						write.close();

					} catch(Exception e) {
						e.printStackTrace();
					}
				}
			} else {
				try {
					Alignment ali= new Alignment(input);
					BufferedWriter write= new BufferedWriter(new FileWriter(new File(output)));
					if (strict) {
						write.write(ali.phylipStrictString() + "\n");
					} else {
						write.write(ali.phylipString() + "\n");
					}
					write.flush();
					write.close();

				} catch(Exception e) {
					e.printStackTrace();
				}
			}


		} else if (operation==3) {
			File inputFiles[] = (new File(inputString)).listFiles();
			for (int i=0;i<inputFiles.length;i++) {
				Alignment ali= new Alignment(CLUSTAL,inputFiles[i]);
				//System.out.println(ali);

				try {
					if (reference!=null) {
						BufferedReader read= new BufferedReader(new FileReader(new File(reference.replace("IN",(inputFiles[i].getName()).substring(0,inputFiles[i].getName().lastIndexOf('.'))))));
						String refString= read.readLine();
						ali.cleanAlignment(refString,inputFiles[i]);
					}
					BufferedWriter write= new BufferedWriter(new FileWriter(new File(output + (inputFiles[i].getName()).substring(0,inputFiles[i].getName().lastIndexOf('.')) + ".phy")));
					write.write(ali.phylipString() + "\n");
					write.flush();
					write.close();

				} catch(Exception e) {
					e.printStackTrace();
				}
			}


		} else if (operation==5) {
			File inputFiles[] = (new File(inputString)).listFiles();
			for (int i=0;i<inputFiles.length;i++) {
				Alignment ali= new Alignment(inputFiles[i]);
				int nbTax= ali.sequences.size();
				int nbSit= ((String)(ali.sequences.elementAt(0))).length();
				System.out.println(inputFiles[i].getName() + "\t" + nbTax + "\t" + nbSit + "\t" + (nbTax*nbTax*nbSit));


			}


		} else if (operation==2) {
			File inputFiles[] = (new File(inputString)).listFiles();
			for (int i=0;i<inputFiles.length;i++) {
				Alignment ali= new Alignment(inputFiles[i]);
				System.out.println(inputFiles[i].getName()+ "\t" + ali.getGapProportion());
			}


		} else if (operation==1) {
			alignment= new Alignment(inputString);
			alignment.cleanAlignment();
			try {
				if (alignment.isValide()) {
					BufferedWriter writer= new BufferedWriter(new FileWriter(new File(output)));
					writer.write(alignment.toString() + "\n");
					writer.flush();
					writer.close();
				}
			} catch(Exception e) {
				e.printStackTrace();
			}
		} else {
			File[] inputFiles = input.listFiles();
			for (int z=0;z<nbBig;z++) {
				int bigger=-1;
				int size=-1;
				for (int i=0;i<inputFiles.length;i++) {
					if (inputFiles[i]!=null) {
						int current=nbLines(inputFiles[i]);
						if (size<current) {
							size=current;
							bigger=i;
						}
					}
				}

				System.out.println(inputFiles[bigger]);
				inputFiles[bigger]=null;
			}
		}

	}
// *******************************************************************************
//STATIC TOOLS
// *******************************************************************************
	public static int nbLines(File f) {
		int res=0;
		try {
			BufferedReader r= new BufferedReader(new FileReader(f));
			String s;
			while ((s=r.readLine())!=null) {
				res++;
			}

			r.close();
		} catch(Exception e) {
			e.printStackTrace();
		}
		//System.out.println(res);
		return res;
	}

// *******************************************************************************
//SUBCLASSES
// *******************************************************************************
	class Alignment {
		Vector names;
		Vector sequences;
		//int nbSequences;
		//int lengthSequences;
		String alignmentName;

		public Alignment(int format,File file) {
					alignmentName= file.getName();
					names= new Vector();
					sequences= new Vector();
					try {
						//read size and length
						BufferedReader read= new BufferedReader(new FileReader(file));


						if (format==CLUSTAL) {

							String s= read.readLine();
							s= read.readLine();
							s= read.readLine();
							s= read.readLine();


							while (s!=null && s.length()>0) {
								StringBuffer nameB= new StringBuffer();
								int i=0;
								while (i<s.length() && s.charAt(i)!=' ') {
									nameB.append(s.charAt(i));
									i++;
								}
								names.addElement(nameB.toString());
								while (i<s.length() && s.charAt(i)==' ') {
									i++;
								}
								StringBuffer buf= new StringBuffer();
								while (i<s.length() && s.charAt(i)!=' ') {
									buf.append(s.charAt(i));
									i++;
								}
								sequences.addElement(buf.toString());
								s= read.readLine();
							}



							s= read.readLine();
							s= read.readLine();
							while (s!=null) {


								int c=0;
								while (s!=null && s.length()>0) {
									//StringBuffer nameB= new StringBuffer();
									int i=0;
									while (i<s.length() && s.charAt(i)!=' ') {
										//nameB.append(s.charAt(i));
										i++;
									}
									//names.addElement(nameB.toString());
									while (i<s.length() && s.charAt(i)==' ') {
										i++;
									}
									StringBuffer buf= new StringBuffer();
									while (i<s.length() && s.charAt(i)!=' ') {
										buf.append(s.charAt(i));
										i++;
									}
									String cur= (String)(sequences.elementAt(c));
									cur=cur+buf.toString();
									sequences.setElementAt(cur,c);
									s= read.readLine();
									c++;

								}




								s= read.readLine();
								s= read.readLine();

							}

						}





						read.close();
						//nbSequences= sequences.size();
						//lengthSequences= ((String)(sequences.elementAt(0))).length();
					} catch(Exception e) {
						e.printStackTrace();

					}


		}

		public Alignment(String input) {
			File file= new File(input);
			names= new Vector();
			sequences= new Vector();
			try {
				//read size and length
				BufferedReader read= new BufferedReader(new FileReader(file));
				String s= read.readLine();
				while (s!=null) {
					names.addElement(s.substring(1,s.length()));
					s= read.readLine();
					StringBuffer buf= new StringBuffer();
					while (s!=null && !s.startsWith(">")) {
						buf.append(s.replace(" ",""));
						s= read.readLine();
					}
					sequences.addElement(new String(buf));
				}

				read.close();
				//nbSequences= sequences.size();
				//lengthSequences= ((String)(sequences.elementAt(0))).length();
			} catch(Exception e) {
				e.printStackTrace();

			}


		}
		public Alignment(File input) {
			names= new Vector();
			sequences= new Vector();
			try {
				//read size and length
				BufferedReader read= new BufferedReader(new FileReader(input));
				String s= read.readLine();
				while (s!=null) {
					names.addElement(s.substring(1,s.length()));
					s= read.readLine();
					StringBuffer buf= new StringBuffer();
					while (s!=null && !s.startsWith(">")) {
						buf.append(s.replace(" ",""));
						s= read.readLine();
					}
					sequences.addElement(new String(buf));
				}

				read.close();
				//nbSequences= sequences.size();
				//lengthSequences= ((String)(sequences.elementAt(0))).length();
			} catch(Exception e) {
				e.printStackTrace();

			}


		}
		public double getGapProportion() {
			double nbChar=0.0;
			double nbBlanks=0.0;
			for (int j=0;j<sequences.size();j++) {
				String sequence= (String)(sequences.elementAt(j));
				for (int i=0;i<sequence.length();i++) {
					nbChar+=1.0;
					char c= sequence.charAt(i);
					if (proteic) {
						if (c=='?' || c =='-') {
							nbBlanks+=1.0;
						}
					} else {
						if (c!='a' && c!='c' && c!='g' && c!='t' && c!='A' && c!='C' && c!='G' && c!='T' && c!='u' && c!='U') {
							nbBlanks+=1.0;
						}

					}
				}

			}
			return nbBlanks/nbChar;
		}
		public void cleanAlignment(String ref, File file) {
			for (int i=0;i<sequences.size();i++) {
				if (ref.indexOf((String)(names.elementAt(i)))==-1) {
					System.out.println("Missing " + (String)(names.elementAt(i)) + " in the " + file.getName() + " file.");
					names.removeElementAt(i);
					sequences.removeElementAt(i);
					i--;
				}
			}

		}

		public void cleanAlignment() {
			boolean cleaned=false;
			while (!cleaned) {
				cleaned=true;
				int res[]= new int[sequences.size()];
				for (int i=0;i<sequences.size();i++) {
					res[i]=0;
					for (int j=0;j<sequences.size();j++) {
						if (i!=j) {
							if (!comparable(i,j)) {
								res[i]++;
								cleaned=false;
							}
						}
					}
				}
				if (!cleaned) {
					int worst=-1;
					int nbWorst=-1;
					for (int i=0;i<res.length;i++) {
						if (nbWorst<res[i]) {
							worst=i;
							nbWorst=res[i];
						}

					}
					names.removeElementAt(worst);
					sequences.removeElementAt(worst);
				}
			}
		}
//java AlignmentEditor -input /auto/dufayard/hogenom/hogenom4/ALN_TREE/ -output /auto/dufayard/hogenom/gamma_beta/ -custom /auto/dufayard/hogenom/gammaproteobacteria_strains.txt /auto/dufayard/hogenom/betaproteobacteria_strains.txt /auto/dufayard/hogenom/stat_gamma_beta.txt

		public boolean filter(Hashtable tab1, Hashtable tab2, BufferedWriter write) {
			int nb1=0;
			int nb2=0;
			int i=0;
			while (i<sequences.size()) {
				String name= (String)(names.elementAt(i));
			/*	if (name.startsWith("HS")) {
					System.out.println(name + tab1.containsKey(name) + tab2.containsKey(name));
				}*/
				String sequence= (String)(sequences.elementAt(i));
				int j=0;
				while (j<name.length() && name.charAt(j)!='0' && name.charAt(j)!='1' && name.charAt(j)!='2' && name.charAt(j)!='3' && name.charAt(j)!='4' && name.charAt(j)!='5' && name.charAt(j)!='6' && name.charAt(j)!='7' && name.charAt(j)!='8' && name.charAt(j)!='9' && name.charAt(j)!='_') {
					j++;
				}
				name= name.substring(0,j);
				if (tab1.containsKey(name)) {
					nb1++;
					StringBuffer buffer= new StringBuffer();
					j=0;
					while (j<sequence.length()) {
						if (sequence.charAt(j)!='-') {
							buffer.append(sequence.charAt(j));
						}
						j++;
					}
					sequences.setElementAt(buffer.toString(),i);
					i++;
				} else if (tab2.containsKey(name)) {
					nb2++;
					StringBuffer buffer= new StringBuffer();
					j=0;
					while (j<sequence.length()) {
						if (sequence.charAt(j)!='-') {
							buffer.append(sequence.charAt(j));
						}
						j++;
					}
					sequences.setElementAt(buffer.toString(),i);
					i++;
				} else {
					names.removeElementAt(i);
					sequences.removeElementAt(i);
				}
			}


			boolean res= ((nb1>=1 && nb2>=3) || (nb2>=1 && nb1>=3));
			if (res) {
				AlignmentEditor.nbValide++;
			}
			if (nb1<3 && nb2<3) {
				AlignmentEditor.nbNo++;
			}
			if (nb1==0 && nb2>2) {
				AlignmentEditor.nbNoGamma++;
			}

			if (nb1>2 && nb2==0) {
				AlignmentEditor.nbNoBeta++;
			}
			try {
				write.write(alignmentName + "\t" + nb1 + "\t" + nb2 + "\n");
				write.flush();
			} catch(Exception e) {
				e.printStackTrace();
			}

			return res;
		}

		public boolean comparable(int a, int b) {
			boolean res=true;
			String seq1= (String)(sequences.elementAt(a));
			String seq2= (String)(sequences.elementAt(b));
			int nbCom=0;
			for (int i=0;i<seq1.length();i++) {
				char c1= seq1.charAt(i);
				char c2= seq2.charAt(i);

				if (c1!='-' && c2!='-' && c1!='?' && c2!='?') {
					nbCom++;
				}

			}
			if (nbCom<seq1.length()/2) {
				res=false;
			}
			return res;
		}

		public boolean isValide() {
			int length= ((String)(sequences.elementAt(0))).length();
			return(length>=2*sequences.size());
		}

		public String phylipStrictString() {
			StringBuffer res= new StringBuffer();
			res.append(sequences.size());
			res.append(" ");
			res.append(((String)(sequences.elementAt(0))).length());
			res.append("\n");
			for (int i=0;i<sequences.size();i++) {
				String name= (String)(names.elementAt(i));
				if (name.length()>9) {
					name= name.substring(0,9);
				}
				res.append(name);
				for (int w=name.length();w<12;w++) {
					res.append(' ');
				}
				res.append((String)(sequences.elementAt(i)));
				res.append("\n");
			}


			return res.toString();
		}

		public String phylipString() {
			StringBuffer res= new StringBuffer();
			res.append(sequences.size());
			res.append(" ");
			res.append(((String)(sequences.elementAt(0))).length());
			res.append("\n");
			for (int i=0;i<sequences.size();i++) {
				res.append((String)(names.elementAt(i)));
				res.append(" ");
				res.append((String)(sequences.elementAt(i)));
				res.append("\n");
			}


			return res.toString();
		}

		public String toString() {
			StringBuffer res= new StringBuffer();
			/*res.append(nbSequences);
			res.append(" sequences of length ");
			res.append(lengthSequences);
			res.append("\n");*/
			for (int i=0;i<sequences.size();i++) {
				res.append(">");
				res.append((String)(names.elementAt(i)));
				res.append("\n");
				res.append((String)(sequences.elementAt(i)));
				res.append("\n");
			}


			return res.toString();
		}







	}


}

