import java.util.*;
import java.io.*;

/**
 * SVG generation script
 * @author Jean-Francois Dufayard
 * @version 1.0
 */
public class ProteinStructureDraw {

// ********************************************************************************************************************
// ***     ATTRIBUTS      ***
// **************************

/**
* Input file
*/
	public static String input;
	
/**
* Output file
*/
	public static String output;	
/**
* parameters
*/
	
	public static int width=800; 
	public static int height=40; 
	public static int margin=20; 
	public static int large=8;
	public static String lineColor="#888888";
	public static int fontSize=14;
	public static String fontFamily="Candara";
	public static int nameMargin=80; 
	public static int legendMargin=70; 
	public static String fontColor="#444444";
	public static int scale=100;
	
	
	public static String[] palette= new String[]{"#3366FF","#33CCCC","#99CC00","#FFCC00","#FF9900","#FF6600","#666699","#969696","#CCFFFF","#CCFFCC","#FFFF99","#99CCFF","#FF99CC","#CC99FF","#FFCC99","#003366","#339966","#003300","#333300","#993300","#993366","#333399","#333333"};



   private static final String NORMAL     = "\u001b[0m";
   private static final String BOLD       = "\u001b[1m";
   private static final String UNDERLINE  = "\u001b[4m";


// ********************************************************************************************************************
// ***     MAIN     ***
// ********************
	public static void main(String[] args) {
		try {
			Hashtable colorTable= new Hashtable();
			Vector colorVector= new Vector();
			String center=null;
			for (int i=0;i<args.length;i=i+2) {
				if (args[i].contains("help")) {
					System.out.println(BOLD);
					System.out.println("NAME:");
					System.out.println(NORMAL);
					System.out.println("\t- Protein domain drawer v1.0 -");
					System.out.println(BOLD);
					System.out.println("SYNOPSIS:");
					System.out.println(NORMAL);
					System.out.println("\tproteindomaindrawer [command args]");
					System.out.println(BOLD);
					System.out.println("MAIN PARAMETERS:");
					System.out.println(BOLD);
					System.out.println("-input" + NORMAL + " " + UNDERLINE + "input_csv_file\n\t" + NORMAL + "The input CSV file containing protein descriptions");
					System.out.println(BOLD);
					System.out.println("-output" + NORMAL + " " + UNDERLINE + "output_svg_file\n\t" + NORMAL + "The output SVG image corresponding to input proteins");
					System.out.println(BOLD);
					System.out.println("-center" + NORMAL + " " + UNDERLINE + "domain_name\n\t" + NORMAL + "Domain that will be aligned between proteins");
					System.out.println(BOLD);
					System.out.println("-color" + NORMAL + " " + UNDERLINE + "domain_name" + NORMAL + " " + UNDERLINE + "color\n\t" + NORMAL + "the domain name or prefix will be colored with the specified color");
					System.out.println(BOLD);
					System.out.println("ADVANCED OPTIONS:");
					System.out.println(BOLD);
					System.out.println("-width" + NORMAL + " " + UNDERLINE + "width\n\t" + NORMAL + "Total width in pixel of the output SVG (default 800)");
					System.out.println(BOLD);
					System.out.println("-height" + NORMAL + " " + UNDERLINE + "height\n\t" + NORMAL + "Height in pixel of one protein display area (default 40)");
					System.out.println(BOLD);
					System.out.println("-margin" + NORMAL + " " + UNDERLINE + "margin\n\t" + NORMAL + "Size in pixel of the white margin around the output SVG (default 20)");
					System.out.println(BOLD);
					System.out.println("-large" + NORMAL + " " + UNDERLINE + "large\n\t" + NORMAL + "Distance in pixel between the line and the border of protein domain rectangles (default 8)");
					System.out.println(BOLD);
					System.out.println("-namemargin" + NORMAL + " " + UNDERLINE + "name_margin\n\t" + NORMAL + "Width in pixel of the specific left margin containing protein names (default 80)");
					System.out.println(BOLD);
					System.out.println("-legendmargin" + NORMAL + " " + UNDERLINE + "legend_margin\n\t" + NORMAL + "Height in pixel of the specific bottom margin containing the legend (default 70)");
					System.out.println(BOLD);
					System.out.println("-scale" + NORMAL + " " + UNDERLINE + "scale\n\t" + NORMAL + "Size in aa or na of the bottom scale (default 100)");
					System.out.println(BOLD);
					System.out.println("-linecolor" + NORMAL + " " + UNDERLINE + "line_color\n\t" + NORMAL + "Color of graphic lines (default 888888)");
					System.out.println(BOLD);
					System.out.println("-fontcolor" + NORMAL + " " + UNDERLINE + "font_color\n\t" + NORMAL + "Color of texts (default 444444)");
					System.out.println(BOLD);
					System.out.println("-fontsize" + NORMAL + " " + UNDERLINE + "font_size\n\t" + NORMAL + "Size of texts in pixel (default 14)");
					System.out.println(BOLD);
					System.out.println("-fontfamily" + NORMAL + " " + UNDERLINE + "font_family\n\t" + NORMAL + "Font family of texts (default Candara)\n\n");
					System.exit(0);
				}
				if (args[i].equalsIgnoreCase("-input")) {
					input= args[i+1];
				}			
				if (args[i].equalsIgnoreCase("-output")) {
					output= args[i+1];
				}			
				if (args[i].equalsIgnoreCase("-center")) {
					center= args[i+1];
				}			
				if (args[i].equalsIgnoreCase("-width")) {
					width= (new Integer(args[i+1])).intValue();;
				}			
				if (args[i].equalsIgnoreCase("-height")) {
					height= (new Integer(args[i+1])).intValue();;
				}				
				if (args[i].equalsIgnoreCase("-margin")) {
					margin= (new Integer(args[i+1])).intValue();;
				}	
				if (args[i].equalsIgnoreCase("-large")) {
					large= (new Integer(args[i+1])).intValue();;
				}			
				if (args[i].equalsIgnoreCase("-fontsize")) {
					fontSize= (new Integer(args[i+1])).intValue();;
				}		
				if (args[i].equalsIgnoreCase("-namemargin")) {
					nameMargin= (new Integer(args[i+1])).intValue();;
				}		
				if (args[i].equalsIgnoreCase("-legendmargin")) {
					legendMargin= (new Integer(args[i+1])).intValue();;
				}			
				if (args[i].equalsIgnoreCase("-scale")) {
					scale= (new Integer(args[i+1])).intValue();;
				}				
				if (args[i].equalsIgnoreCase("-linecolor")) {
					lineColor= args[i+1].replace("__pd__","#");
				}			
				if (args[i].equalsIgnoreCase("-fontfamily")) {
					fontFamily= args[i+1];
				}			
				if (args[i].equalsIgnoreCase("-fontcolor")) {
					fontColor= args[i+1].replace("__pd__","#");
				}						
				if (args[i].equalsIgnoreCase("-color")) {
					colorTable.put(args[i+1],args[i+2].replace("__pd__","#"));
					colorVector.addElement(args[i+1]);
					i++;
				}
			}
			int max=0;
			BufferedReader read= new BufferedReader(new FileReader(new File(input)));
			BufferedWriter write= new BufferedWriter(new FileWriter(new File(output)));
			
			Hashtable proteinTable= new Hashtable();
			Vector proteinVector= new Vector();
			Hashtable centers= new Hashtable();
			
			// parsing CSV input file
			String s = read.readLine();
			s = read.readLine();
			double maxCenter=0;
			while (s!=null) {
				String[] splits=s.split(";");
				Vector thisVector=null;
				if (proteinTable.containsKey(splits[0])) {
					thisVector=(Vector)(proteinTable.get(splits[0]));
				} else {
					thisVector=new Vector();
					proteinTable.put(splits[0],thisVector);
					proteinVector.addElement(splits[0]);
				}
				if (splits[1].equalsIgnoreCase(center)) {
					double localCenter=(new Double(((new Double(splits[3])).doubleValue()+(new Double(splits[2])).doubleValue())/2.0)).doubleValue();
					if (localCenter>maxCenter) {
						maxCenter=localCenter;
					}
					//System.out.println(new Double(localCenter));
					centers.put(splits[0],localCenter);
				} else {
					if (!centers.containsKey(splits[0])) {
						centers.put(splits[0],0.0);			
					}
				}
				if (splits[1].length()==0) {
					splits[1]="globallength";
				}
				thisVector.addElement(splits[1]);
				thisVector.addElement(splits[2]);
				thisVector.addElement(splits[3]);
				s = read.readLine();
			}
			read.close();
				
			
			write.write("<?xml version=\"1.0\" standalone=\"yes\"?>\n<svg version=\"1.1\" width=\"");
			write.write((new Integer(width)).toString());				
			write.write("\" height=\"");
			write.write((new Integer(height*proteinVector.size()+2*margin+legendMargin)).toString());
			write.write("\" fill=\"none\" stroke=\"none\" stroke-linecap=\"square\" stroke-miterlimit=\"10\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n");
			write.flush();
			
			//compute max coordinate
			double maxEnd=0.0;
			for (int i=0;i<proteinVector.size();i++) {
				Vector thisVector=(Vector)(proteinTable.get((String)(proteinVector.elementAt(i))));
				double localCenter=((Double)(centers.get((String)(proteinVector.elementAt(i))))).doubleValue();
				int j=0;
				try {
					while (!((String)(thisVector.elementAt(j))).equals("globallength")) {
						j++;
					}
				} catch(Exception exp) {
					System.out.println("Missing overall length of " + (String)(proteinVector.elementAt(i)));
					System.exit(1);
				}
				double localEnd=0.0;
				if (localCenter==0.0) {
					localEnd=(new Double((String)(thisVector.elementAt(j+2)))).doubleValue();
				} else { 
					localEnd=(new Double((String)(thisVector.elementAt(j+2)))).doubleValue() + maxCenter-localCenter;
				}
				if (localEnd>maxEnd) {
					maxEnd=localEnd;
				}
				//System.out.println((String)(thisVector.elementAt(j+1) + ":" + (String)(thisVector.elementAt(j+2))));
			
			
			}
			double pToAA= ((double)width-2*(double)margin-nameMargin)/maxEnd;
			write.write("<line x1=\"" + (int)((double)margin+nameMargin+maxCenter*pToAA) + "\" y1=\"" + margin + "\" x2=\"" + (int)((double)margin+nameMargin+maxCenter*pToAA) + "\" y2=\"" + (height*proteinVector.size()+margin) + "\" style=\"stroke:" + lineColor + ";stroke-width:1;stroke-dasharray: 3, 3;\"/>\n");
			
			int colorJogger=0;	
			for (int i=0;i<proteinVector.size();i++) {
				String thisName=(String)(proteinVector.elementAt(i));
				Vector thisVector=(Vector)(proteinTable.get(thisName));
				double localCenter=((Double)(centers.get(thisName))).doubleValue();
				int j=0;
				while (!((String)(thisVector.elementAt(j))).equals("globallength")) {
					j++;
				}
				double delta=0.0;
				if (localCenter!=0.0) {
					delta=maxCenter-localCenter;
				}
				double end=(new Double((String)(thisVector.elementAt(j+2)))).doubleValue();
				end= (double)margin+(double)nameMargin+(delta+end)*pToAA;
				double start=(double)margin+(double)nameMargin+(delta)*pToAA;
				//System.out.println(start + "::" + end);
				write.write("<line x1=\"" + (int)start + "\" y1=\"" + (margin+height*i+height/2) + "\" x2=\"" + (int)end + "\" y2=\"" + (margin+height*i+height/2) + "\" style=\"stroke:" + lineColor + ";stroke-width:1;\"/>\n");
				write.write("<text x=\"" + (int)margin + "\" y=\"" + (margin+height*i+height/2+fontSize/2-2) + "\" style=\"font-family:" + fontFamily + ";font-size:" + fontSize + ";fill:" + fontColor + ";\">" + thisName + "</text>\n");
				write.flush();
				j=0;
				while (j<thisVector.size()) {
					if (!((String)(thisVector.elementAt(j))).equals("globallength")) {
						String domName=(String)(thisVector.elementAt(j));
						String color="#CCCCCC";
						int k=0;
						while (k<colorVector.size() && domName.indexOf((String)(colorVector.elementAt(k)))==-1) {
							k++;
						}
						if (k<colorVector.size()) {
							color=(String)(colorTable.get((String)(colorVector.elementAt(k))));
						} else {
							if (colorJogger>=palette.length) {
								color="#000000";
							} else {
								color=palette[colorJogger];
								
							}
							colorJogger++;
							if (!colorTable.containsKey(domName)) {
								colorVector.addElement(domName);
								colorTable.put(domName,color);
							}
						
						}
					
						start=(new Double((String)(thisVector.elementAt(j+1)))).doubleValue();
						start= (double)margin+(double)nameMargin+(delta+start)*pToAA;
						end=(new Double((String)(thisVector.elementAt(j+2)))).doubleValue();
						end= (double)margin+(double)nameMargin+(delta+end)*pToAA;
						if ((int)(end-start)<1) {		 			
							write.write("<rect x=\"" + (int)start + "\" y=\"" + (margin+height*i+height/2-large) + "\" width=\"" + 1 + "\" height=\"" + (large*2+1) + "\" style=\"fill:" + color + ";stroke:rgb(0,0,0);stroke-width:0;\"/>\n");
						} else { 			
							write.write("<rect x=\"" + (int)start + "\" y=\"" + (margin+height*i+height/2-large) + "\" width=\"" + (int)(end-start) + "\" height=\"" + (large*2+1) + "\" style=\"fill:" + color + ";stroke:rgb(0,0,0);stroke-width:0;\"/>\n");
						}
						write.flush();	
					}			
				
				
					j+=3;
				}				
			
			}
			
			Collections.sort((List)colorVector);
			// display legend
			int deltaX=0;
			for (int i=0;i<colorVector.size();i++) {
				String thisName= (String)(colorVector.elementAt(i));
				String thisColor= (String)(colorTable.get(thisName));
				write.write("<rect x=\"" + (margin+deltaX) + "\" y=\"" + (height*proteinVector.size()+margin+legendMargin-large) + "\" width=\"" + (4*large) + "\" height=\"" + (large*2+1) + "\" style=\"fill:" + thisColor + ";stroke:rgb(0,0,0);stroke-width:0;\"/>\n");
				write.write("<text x=\"" + (margin+deltaX+4*large+5) + "\" y=\"" + (height*proteinVector.size()+margin+legendMargin+fontSize/2-2) + "\" style=\"font-family:" + fontFamily + ";font-size:" + fontSize + ";fill:" + fontColor + ";\">" + thisName + "</text>\n");
				deltaX+=4*large + thisName.length()*fontSize/3+40;
			}

			//display scale
			write.write("<text text-anchor=\"middle\" x=\"" + (margin+nameMargin) + "\" y=\"" + (height*proteinVector.size()+margin+large+8+fontSize/2+2+fontSize) + "\" style=\"font-family:" + fontFamily + ";font-size:" + fontSize + ";fill:" + fontColor + ";\">" + 0 + "</text>\n");
			write.write("<line x1=\"" + (margin+nameMargin) + "\" y1=\"" + (height*proteinVector.size()+margin+large+fontSize) + "\" x2=\"" + (margin+nameMargin) + "\" y2=\"" + (height*proteinVector.size()+margin+large+5+fontSize) + "\" style=\"stroke:" + lineColor + ";stroke-width:1;\"/>\n");			
			int i=scale;
			while (i<maxEnd) {
				write.write("<text text-anchor=\"middle\" x=\"" + (margin+nameMargin+i*pToAA) + "\" y=\"" + (height*proteinVector.size()+margin+large+8+fontSize/2+2+fontSize) + "\" style=\"font-family:" + fontFamily + ";font-size:" + fontSize + ";fill:" + fontColor + ";\">" + i + "</text>\n");
				write.write("<line x1=\"" + (margin+nameMargin+i*pToAA) + "\" y1=\"" + (height*proteinVector.size()+margin+large+fontSize) + "\" x2=\"" + (margin+nameMargin+i*pToAA) + "\" y2=\"" + (height*proteinVector.size()+margin+large+5+fontSize) + "\" style=\"stroke:" + lineColor + ";stroke-width:1;\"/>\n");			
				i+=scale;
			}
			write.write("<line x1=\"" + (margin+nameMargin) + "\" y1=\"" + (height*proteinVector.size()+margin+large+fontSize) + "\" x2=\"" + (width-margin) + "\" y2=\"" + (height*proteinVector.size()+margin+large+fontSize) + "\" style=\"stroke:" + lineColor + ";stroke-width:1;\"/>\n");	
			write.write("<text text-anchor=\"middle\" x=\"" + (width-margin) + "\" y=\"" + (height*proteinVector.size()+margin+large+8-fontSize/2) + "\" style=\"font-family:" + fontFamily + ";font-size:" + fontSize + ";fill:" + fontColor + ";\">" + ((int)(maxEnd)) + "</text>\n");
			write.write("<line x1=\"" + (width-margin) + "\" y1=\"" + (height*proteinVector.size()+margin+large+fontSize) + "\" x2=\"" + (width-margin) + "\" y2=\"" + (height*proteinVector.size()+margin+large-5+fontSize) + "\" style=\"stroke:" + lineColor + ";stroke-width:1;\"/>\n");			
			write.write("</svg>");
			write.flush();
			
			write.close();
		

		} catch(Exception e) {
			e.printStackTrace();
		}
	}







}