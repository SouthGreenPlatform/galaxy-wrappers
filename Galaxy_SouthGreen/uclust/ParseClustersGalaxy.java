import java.util.*;
import java.io.*;

/**
 * NAME:
 * <p>
 * 		- Parse clusters -
 * <p>
 * SYNOPSIS:
 * <p>
 * 		parseclusters [command args]
 * <p>
 * OPTIONS:
 * <p>
 * 		-sequences sequences_file :
 *
 * 			File containing every sequences in fasta format
 * <p>
 * 		-clusters clusters_file :
 *
 * 			File containing every cluster, one per line, coma or space separated
 * <p>
 * 		-output output_directory :
 *
 * 			The directory to stock output fasta files
 *
 * @author Jean-Francois Dufayard
 * @version 1.0
 */
public class ParseClustersGalaxy {

// ********************************************************************************************************************
// ***     ATTRIBUTS      ***
// **************************

/**
* Sequences file
*/
	public static File sequencesFile;

/**
* Clusters file
*/
	public static File clustersFile;

/**
* Output statistic file
*/
	public static File statisticFile;

/**
* Output file
*/
	public static String output;

	private static final String NORMAL     = "\u001b[0m";
	private static final String BOLD       = "\u001b[1m";
	private static final String UNDERLINE  = "\u001b[4m";

/**
*
* /usr/local/jdk1.6.0_20/bin/java -mx4096m -cp /home/dufayard/utilities/ ParseClusters -sequences sativa_sorted.fasta -clusters results_0.5.uc -stats stats_0.5_sativaonly.csv -trans -output /home/dufayard/uclust/clusters_0.5_sativaonly/
*
* /usr/local/jdk1.6.0_20/bin/java -mx4096m -cp /home/dufayard/utilities/ ParseClusters -sequences sativa_sorted.fasta -clusters results_0.4.uc -stats stats_0.4.csv -trans -output /home/dufayard/uclust/clusters_0.4/ -complement RC1.fasta RC1 RC2.fasta RC2 RC3.fasta RC3 RC4.fasta RC4 RC5.fasta RC5 RC6.fasta RC6 RC7.fasta RC7 RC8.fasta RC8 RC9.fasta RC9 RC10.fasta RC10 RS7.fasta RS7
*
* /usr/local/jdk1.6.0_20/bin/java -mx4096m -cp /home/dufayard/utilities/ ParseClusters -sequences sativa_sorted.fasta -clusters results_0.5.uc -stats stats_0.5.csv -trans -output /home/dufayard/uclust/clusters_0.5/ -complement RC1.fasta RC1 RC2.fasta RC2 RC3.fasta RC3 RC4.fasta RC4 RC5.fasta RC5 RC6.fasta RC6 RC7.fasta RC7 RC8.fasta RC8 RC9.fasta RC9 RC10.fasta RC10 RS7.fasta RS7
*
* /usr/local/jdk1.6.0_20/bin/java -mx4096m -cp /home/dufayard/utilities/ ParseClusters -sequences sativa_sorted.fasta -clusters results_0.6.uc -stats stats_0.6.csv -trans -output /home/dufayard/uclust/clusters_0.6/ -complement RC1.fasta RC1 RC2.fasta RC2 RC3.fasta RC3 RC4.fasta RC4 RC5.fasta RC5 RC6.fasta RC6 RC7.fasta RC7 RC8.fasta RC8 RC9.fasta RC9 RC10.fasta RC10 RS7.fasta RS7
*
* /usr/local/jdk1.6.0_20/bin/java -mx4096m -cp /home/dufayard/utilities/ ParseClusters -sequences sativa_sorted.fasta -clusters results_0.7.uc -stats stats_0.7.csv -trans -output /home/dufayard/uclust/clusters_0.7/ -complement RC1.fasta RC1 RC2.fasta RC2 RC3.fasta RC3 RC4.fasta RC4 RC5.fasta RC5 RC6.fasta RC6 RC7.fasta RC7 RC8.fasta RC8 RC9.fasta RC9 RC10.fasta RC10 RS7.fasta RS7
*/

// ********************************************************************************************************************
// ***     MAIN     ***
// ********************
	public static void main(String[] args) {
		try {
			boolean verbose=false;
			boolean trans=false;

			Vector complementFiles= new Vector();
			Vector complementTaxa= new Vector();

			for (int i=0;i<args.length;i=i+2) {
				if (args[i].contains("help")) {
					System.out.println(BOLD);
					System.out.println("NAME:");
					System.out.println(NORMAL);
					System.out.println("\t- Parse clusters -");
					System.out.println(BOLD);
					System.out.println("SYNOPSIS:");
					System.out.println(NORMAL);
					System.out.println("\tparseclusters [command args]");
					System.out.println(BOLD);
					System.out.println("OPTIONS:");
					System.out.println(BOLD);
					System.out.println("-sequences" + NORMAL + " " + UNDERLINE + "sequences_file\n\t" + NORMAL + "File containing every sequences in fasta format");
					System.out.println(BOLD);
					System.out.println("-clusters" + NORMAL + " "  + UNDERLINE + "clusters_file\n\t" + NORMAL + "File containing every cluster, one per line, coma or space separated ; or standard '.uc' UClust output file");
					System.out.println(BOLD);
					System.out.println("-stats" + NORMAL + " "  + UNDERLINE + "output_statistic_file\n\t" + NORMAL + "File to stock statistics");
					System.out.println(BOLD);
					System.out.println("-trans\n\t" + NORMAL + "Add transcript numbers to the stats file");
					System.out.println(BOLD);
					System.out.println("-output" + NORMAL + " "  + UNDERLINE + "output_directory\n\t" + NORMAL + "The directory to stock output fasta files");
					System.exit(0);
				}

				if (args[i].equalsIgnoreCase("-sequences")) {
					sequencesFile= new File(args[i+1]);
				}
				if (args[i].equalsIgnoreCase("-clusters")) {
					clustersFile= new File(args[i+1]);
				}
				if (args[i].equalsIgnoreCase("-output")) {
					output= args[i+1];
				}
				if (args[i].equalsIgnoreCase("-verbose")) {
					i--;
					verbose=true;
				}
			}
			System.out.print("Parsing sequences...");

			BufferedReader read= new BufferedReader(new FileReader(sequencesFile));
			String s= read.readLine();
			Hashtable sequences= new Hashtable();
			int nbSequences=0;
			while (s!=null) {
				String key= s.substring(1,s.length());
				s= read.readLine();
				StringBuffer seq= new StringBuffer();
				while (s!=null && !s.startsWith(">")) {
					seq.append(s);
					s= read.readLine();
				}
				if (!sequences.containsKey(key)) {
					sequences.put(key,seq.toString());
					nbSequences++;
				}


			}
			System.out.println(" Done.");
			System.out.println("Number of sequences: " + nbSequences);

			read.close();


			Hashtable complements= new Hashtable();


			read= new BufferedReader(new FileReader(clustersFile));
			Hashtable clusters = new Hashtable();
			Vector clusterIds= new Vector();
			s= read.readLine();
			int nbClusters=0;
			int nbClassifiedSequences=0;
			Hashtable trace = new Hashtable();
			while (s!=null) {
				String words[]= s.split("\t");

				if (words[0].equals("S")) {
					if (trace.containsKey(words[words.length-2])) {
						if (verbose) {
							String oldId= (String)(trace.get(words[words.length-2]));
							System.out.println("Wants to insert " + words[words.length-2] + " in cluster " + words[1] + ", but is already inserted in cluster " + oldId + ".");
						}
					} else {
												Vector newCluster= new Vector();
						newCluster.addElement(words[words.length-2]);
						clusterIds.addElement(words[1]);
						clusters.put(words[1],newCluster);
						trace.put(words[words.length-2],words[1]);
						//System.out.println("ID=" + words[1] + " SEQ=" + words[words.length-2]);
					}

				} else if (words[0].equals("H")) {
					Vector currentCluster= (Vector)(clusters.get(words[1]));
					for (int i=8;i<words.length;i++) {
						if (trace.containsKey(words[i])) {
							if (verbose) {
								String oldId= (String)(trace.get(words[i]));
								System.out.println("Wants to insert " + words[i] + " in cluster " + words[1] + ", but is already inserted in cluster " + oldId + ".");
							}
						} else {
							trace.put(words[i],words[1]);
							//System.out.println("ID=" + words[1] + " SEQ=" + words[i]);
							if (currentCluster.size()==1) {
								nbClassifiedSequences+=2;
								nbClusters++;
							} else {
								nbClassifiedSequences++;
							}
							currentCluster.addElement(words[i]);
						}
					}
				}
				s= read.readLine();
			}
			System.out.println("Number of clusters: " + nbClusters + "\nNumber of clustered sequences: " + nbClassifiedSequences);

			BufferedWriter write= new BufferedWriter(new FileWriter(output));

			for (int i=0;i<clusterIds.size();i++) {
				String currentId=(String)(clusterIds.elementAt(i));

				Vector current= (Vector)(clusters.get(currentId));


				write.write("\nCluster " + (i+1) + ":\n\n");
				write.flush();

				int clusterSize=0;

				for (int j=0;j<current.size();j++) {
					String currentSequence=(String)(current.elementAt(j));
					write.write(">" + currentSequence + "\n");
					write.flush();
					write.write(((String)(sequences.get(currentSequence))) + "\n");
					write.flush();
					clusterSize++;


					if (complementFiles.size()>0) {
						//complete the file with complementary taxa
						for (int w=0;w<complementTaxa.size();w++) {
							String taxon= (String)(complementTaxa.elementAt(w));
							Hashtable complement= (Hashtable)(complements.get(taxon));
							String currentKey= currentSequence.substring(0,currentSequence.lastIndexOf("_")) + "_"  + taxon;
							if (complement.containsKey(currentKey)) {
								write.write(">" + currentKey + "\n");
								write.flush();
								write.write(((String)(complement.get(currentKey))) + "\n");
								write.flush();
								clusterSize++;
							}

						}




					}




				}



			}
			write.close();


			read.close();
		} catch(Exception e) {
			e.printStackTrace();
		}
	}







}