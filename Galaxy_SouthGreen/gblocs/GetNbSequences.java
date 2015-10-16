import java.util.*;
import java.io.*;

public class GetNbSequences {

// ********************************************************************************************************************
// ***     ATTRIBUTS      ***
// **************************

/**
* Sequences file
*/
	public static File input;




// ********************************************************************************************************************
// ***     MAIN     ***
// ********************
	public static void main(String[] args) {
		try {
			input= new File(args[0]);

			BufferedReader read= new BufferedReader(new FileReader(input));
			String s= read.readLine();
			Hashtable sequences= new Hashtable();
			Vector ids= new Vector();
			boolean fasta=true;
			if (s.startsWith(">")) {
				// Fasta mod
				while (s!=null) {
					String key= s.substring(1,s.length());
					s= read.readLine();
					StringBuffer seq= new StringBuffer();
					while (s!=null && !s.startsWith(">")) {
						seq.append(s);
						s= read.readLine();
					}
					sequences.put(key,seq.toString());
					ids.addElement(key);

				}

			} else {
				fasta=false;
				// Phylip mod
				StringBuffer nbBuffer= new StringBuffer();
				StringBuffer lenBuffer= new StringBuffer();
				int i=0;
				while (s.charAt(i)==' ' || s.charAt(i)=='\t') {
					i++;
				}
				while (s.charAt(i)!=' ' && s.charAt(i)!='\t') {
					nbBuffer.append(s.charAt(i));
					i++;
				}
				while (s.charAt(i)==' ' || s.charAt(i)=='\t') {
					i++;
				}
				while (i<s.length() && s.charAt(i)!=' ' && s.charAt(i)!='\t') {
					lenBuffer.append(s.charAt(i));
					i++;
				}
				s= read.readLine();
				boolean first=true;
				while (s!=null) {
					int k=0;
					while (s!=null && s.length()>2) {
						i=0;
						StringBuffer key= new StringBuffer();
						StringBuffer seq= new StringBuffer();
						while (s.charAt(i)==' ' || s.charAt(i)=='\t') {
							i++;
						}
						if (first) {
							while (s.charAt(i)!=' ' && s.charAt(i)!='\t') {
								key.append(s.charAt(i));
								i++;
							}
							while (s.charAt(i)==' ' || s.charAt(i)=='\t') {
								i++;
							}
						}
						while (i<s.length()) {
							if (s.charAt(i)!=' ') {
								seq.append(s.charAt(i));
							}
							i++;
						}
						if (first) {
							ids.addElement(key.toString());
							sequences.put(key.toString(),seq.toString());
						} else {
							String start= (String)(sequences.get((String)(ids.elementAt(k))));
							sequences.put(((String)(ids.elementAt(k))),start+seq.toString());
						}
						s= read.readLine();
						k++;
					}
					while (s!=null && s.length()<2) {
						s= read.readLine();
					}
					first=false;
				}
			}

			read.close();
			System.out.println(ids.size());

		} catch(Exception e) {
			e.printStackTrace();
		}
	}







}