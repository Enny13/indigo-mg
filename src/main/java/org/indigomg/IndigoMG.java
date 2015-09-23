package org.indigomg;
import java.io.IOException;
import java.util.*;
import com.ggasoftware.indigo.Indigo;
import com.ggasoftware.indigo.IndigoObject;

/**
 * Molecule Generation
 * The main class collecting parameters and setting global objects
 *
 * @author Artem Malykh
 */
public class IndigoMG {
	private static String fragments = null;
	private static Indigo indigo = new Indigo();

	public static void main(String[] args) throws IOException{
		String formula = null;
		String out = "default_out.sdf";

		if (args.length > 0) {

			for(int i = 0; i < args.length; i++){
				if(args[i].equals("-h")){
					System.out.println("indigo-mg generates chemical structures");
					System.out.println("");
					System.out.println("Usage: java -jar indigo-mg.jar -ec <elemental_composition> [-o <out_file.sdf>, -fr <in_fragments.sdf>]");
					System.out.println("");
					System.out.println("Required Parameters");
					System.out.println("-ec:  elemental composition of the molecules to be generated.");
					System.out.println("");
					System.out.println("Optional Parameters");
					System.out.println("-o:   SDF file where to store the molecules. ");
					System.out.println("-fr:  SDF file containing prescribed one or multiple substructures. In the case");
					System.out.println("         of multiple substructures, they have to be non-overlapping. ");
					System.out.println("");
					System.out.println("");
					System.out.println("Examples:");
					System.out.println("java -jar indigo-mg.jar -ec C6H6");
					System.out.println("");
					System.out.println("java -jar indigo-mg.jar -ec C6H6 -o out_C6H6.sdf");
					System.out.println("");
					System.out.println("java -jar indigo-mg.jar -ec C2H5NO2 -fr fragment_CO2.sdf");
					System.out.println("");

					System.exit(1);
				}
				if(args[i].equals("-ec")){
					try {
						formula = args[i + 1];
					} catch (Exception e) {
						System.err.println("No formula provided");
						System.exit(1);
					}
				}
				else if(args[i].equals("-o")){
					try {
						out = args[i + 1];
					} catch (Exception e) {
						System.err.println("No output file provided");
						System.exit(1);
					}
				}
				else if(args[i].equals("-fr")){
					try {
						fragments = args[i + 1];
					} catch (Exception e) {
						System.err.println("No file with prescribed substructures provided");
						System.exit(1);
					}
				}
			}
		}
		else{
			System.err.println("Provide at least an elemental composition. Type indigo-mg.jar -h for help");
			System.exit(1);
		}


		long before = System.currentTimeMillis();

		List<IndigoObject> res = null;
		try {
			res = new MoleculesGenerator(indigo).generateMolecules(formula, fragments);
			IndigoObject saver = indigo.writeFile(out);
			for (IndigoObject mol : res) {
				saver.sdfAppend(mol);
			}
		} catch (FragmentsOverlapException e) {
		}
		System.out.println("Total mols: "+ res.size());
		System.out.println("Took " + (System.currentTimeMillis() - before) + " ms.");

	}
}
