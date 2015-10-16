package rapgreen;
import java.util.*;
import java.io.*;

/**
 * @author Jean-Francois Dufayard
 * @version 1.0
 * Filter PhyloXML trees, removing specified taxa
 */
public class TreeFilter {

// ********************************************************************************************************************
// ***     ATTRIBUTS      ***
// **************************
/**
* Input tree file
*/
	public static File treeFile;
	
/**
* Input species tree file
*/
	public static File speciesFile;
	
/**
* Output cluster file
*/
	public static File outputFile;
	
/**
* List of taxids to remove
*/
	public static Vector idList;
	
/**
* Table of taxids to remove
*/
	public static Hashtable idTable;
	
/**
* Support threshold
*/

	//public static double support=null;
	
	public static boolean clean=false;
	
	public static boolean newick=false;
	
	public static boolean newick2NHX=false;
	
	public static String key="";

// ********************************************************************************************************************
// ***     MAIN     ***
// ********************
	public static void main(String[] args) {
		outputFile=null;
		treeFile=null;
		idList= new Vector();
		idTable= new Hashtable();
		try {
			for (int i=0;i<args.length;i=i+2) {
				if (args[i].equalsIgnoreCase("-input")) {
					treeFile= new File(args[i+1]);
				}
				if (args[i].equalsIgnoreCase("-output")) {
					outputFile= new File(args[i+1]);
				}
				if (args[i].equalsIgnoreCase("-remove")) {
					idList.addElement(args[i+1]);
					idTable.put(args[i+1],"");
				}
				if (args[i].equalsIgnoreCase("-clean")) {
					clean=true;
					i--;
				}
				if (args[i].equalsIgnoreCase("-newick")) {
					newick=true;
					i--;
				}
				if (args[i].equalsIgnoreCase("-key")) {
					key= args[i+1];
				}
				if (args[i].equalsIgnoreCase("-newick2NHX")) {
					newick2NHX=true;
					speciesFile= new File(args[i+1]);
					
				}
								
			}
			
			if (newick2NHX) {
 				File[] treeFiles = treeFile.listFiles();	

				TreeReader read= null;
				BufferedReader buf= new BufferedReader(new FileReader(speciesFile));
				String test= buf.readLine();
				buf.close();
				if (test.endsWith(";")) {
					read= new TreeReader(speciesFile,TreeReader.NEWICK);
				} else {
					read= new TreeReader(speciesFile,TreeReader.XML);
				}
				Tree speciesTree= read.nextTree();
				speciesTree.pretreatment();
			
	        	for (int i=0;i<treeFiles.length;i++) {		
	        		if (treeFiles[i].getName().contains(key)) {
						read= null;
						buf= new BufferedReader(new FileReader(treeFiles[i]));
						test= buf.readLine();
						buf.close();
						if (test.endsWith(";")) {
							read= new TreeReader(treeFiles[i],TreeReader.NEWICK);
						} else {
							read= new TreeReader(treeFiles[i],TreeReader.XML);
						}
						Tree tree= read.nextTree();
						tree.pretreatment();
						BufferedWriter write= new BufferedWriter(new FileWriter(new File(outputFile.getPath() + "/" + treeFiles[i].getName() + ".nhx")));
						write.write(tree.getNHXNewick(speciesTree) + "\n");	
						write.flush();
						write.close();							
						
					}	        	
	        	
	        	
	        	
	        	}		
			
				System.exit(0);
			
			} else if (clean) {
				TreeReader reader= new TreeReader(treeFile,TreeReader.NEWICK);
				Tree tree= reader.nextTree();	
				tree.pretreatment();
				
				for (int i=0;i<tree.leafVector.size();i++) {
					Tree leaf= (Tree)(tree.leafVector.elementAt(i));	
					leaf.label= leaf.label.substring(0,leaf.label.indexOf("_"));
					while (leaf.label.startsWith("0")) {
						leaf.label= leaf.label.substring(1,leaf.label.length());
						
						
					}
					
				}
				String phylo= tree.getNewick();
				if (outputFile!=null) {
					BufferedWriter write= new BufferedWriter(new FileWriter(outputFile));
					write.write(phylo+"\n");
					write.flush();
					write.close();
				} else {
					System.out.println(phylo);				
				}
				System.exit(0);	
			}
			//System.out.println(tree);
			
			TreeReader reader= new TreeReader(treeFile,TreeReader.XML);
			Tree tree= reader.nextTree();	
			Hashtable excluded= new Hashtable();
			tree.taxonomicPretreatment();
			//System.out.println(tree.getNewick());
			/*if (support!=null) {
				tree.collapseSupport(support);
				String phylo= tree.getNewick();
				if (outputFile!=null) {
					BufferedWriter write= new BufferedWriter(new FileWriter(outputFile));
					write.write(phylo+"\n");
					write.flush();
					write.close();
				} else {
					System.out.println(phylo);				
				}
				
			} else {*/

			//}

		} catch(Exception e) {
			e.printStackTrace();
			System.out.println("Usage for standard tree filtering:\njava -jar TreeFilter.jar -input your_tree_file -remove taxid_to_remove [-remove taxid_to_remove] [-output your_output_file] [-newick] [-newick2NHX] [-key]\n");
		}
	}

// ********************************************************************************************************************
// ***     PRIVATE CONSTRUCTOR TOOLS     ***
// *****************************************
	private static void filterRelations(StringBuffer relations, Hashtable excluded) {
		String s=null;
		try {
			BufferedReader read= new BufferedReader(new FileReader(treeFile));
			s= read.readLine();
			while (s!=null && !s.contains("<sequence_relation")) {
				s= read.readLine();
			}
			while (s!=null) {
				StringBuffer id1= new StringBuffer();
				StringBuffer id2= new StringBuffer();
				int i=0;
				while (s.charAt(i)!='"') {
					i++;				
				}
				i++;
				while (s.charAt(i)!='"') {
					id1.append(s.charAt(i));
					i++;				
				}
				i++;	
				while (s.charAt(i)!='"') {
					i++;				
				}
				i++;
				while (s.charAt(i)!='"') {
					id2.append(s.charAt(i));
					i++;				
				}
				if (excluded.containsKey(id1.toString()) || excluded.containsKey(id2.toString())) {
					s=read.readLine();					
					while (s!=null && !s.contains("<sequence_relation")) {
						s= read.readLine();
					}
				} else {
					relations.append(s + "\n");
					s=read.readLine();					
					while (s!=null && !s.contains("</sequence_relation")) {						
						relations.append(s + "\n");
						s= read.readLine();
					}						
					relations.append(s + "\n");
					while (s!=null && !s.contains("<sequence_relation")) {
						s= read.readLine();
					}
					
				}
				
				
			}
			
			
			
								
			read.close();
		} catch(Exception e) {
			System.out.println("bugged line: " + s);
			e.printStackTrace();
			
		}
	}



}