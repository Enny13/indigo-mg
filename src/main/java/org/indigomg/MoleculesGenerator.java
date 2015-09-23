package org.indigomg;

import com.ggasoftware.indigo.Indigo;
import com.ggasoftware.indigo.IndigoObject;
import org.apache.log4j.Logger;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class used for generating list of molecules with given elemental composition
 * Created by Artem Malykh on 16.09.15.
 */
public class MoleculesGenerator {

    public final static Map<String, Integer> VALENCES;
    private static Logger logger = Logger.getLogger(MoleculesGenerator.class);

    private static final int MAX_OPENINGS_INDEX = 0;
    private static final int NH_INDEX = 0;

    static {
        VALENCES = new HashMap<>();
        VALENCES.put("C", 4);
        VALENCES.put("N", 3);
        VALENCES.put("O", 2);
        VALENCES.put("S", 6);
        VALENCES.put("P", 5);
        VALENCES.put("F", 1);
        VALENCES.put("I", 1);
        VALENCES.put("Cl", 1);
        VALENCES.put("Br", 1);
    }

    /**
     * Encapsulates source and destination indexes
     */
    public static class SD {
        int s;
        int d;

        public SD(int s, int d) {
            this.s = s;
            this.d = d;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            SD sd = (SD) o;

            return s == sd.s && d == sd.d;

        }

        @Override
        public int hashCode() {
            int result = s;
            result = 31 * result + d;
            return result;
        }
    }

    private Indigo indigo;

    public MoleculesGenerator(Indigo indigo) {
        this.indigo = indigo;
    }

    /**
     * Generate all possible molecules with given elemental composition and containing given fragments
     * @param elementalComposition Elemental composition
     * @param fragments Fragments separated by space
     * @return List of possible molecules
     */
    public List<IndigoObject> generateMolecules(String elementalComposition, List<IndigoObject> fragments) throws FragmentsOverlapException {
        IndigoObject container = loadElementalComposition(elementalComposition);
        return generateMolecules(container, fragments);
    }


    /**
     * Generate all possible molecules with given elemental composition and containing given fragments
     * @param elementalComposition Elemental composition
     * @param fragments Fragments separated by space
     * @return List of possible molecules
     */
    public List<IndigoObject> generateMolecules(String elementalComposition, String fragments) throws FragmentsOverlapException {
        IndigoObject container = loadElementalComposition(elementalComposition);


        List<IndigoObject> fs = null;
        if (fragments != null) {
            String[] frags = fragments.split(" ");
            fs = new LinkedList<>();
            for (String frag : frags) {
                fs.add(indigo.loadMolecule(frag));
            }
        }
        return generateMolecules(container, fs);
    }

    /**
     * Generate all possible molecules with given elemental composition and containing given fragments.
     * @param container Elemental composition
     * @param fragments List of non-overlapping fragments
     * @return List of possible molecules
     * @throws FragmentsOverlapException
     */
    public List<IndigoObject> generateMolecules(IndigoObject container, List<IndigoObject> fragments) throws FragmentsOverlapException {
        List<IndigoObject> res = new LinkedList<>();
        if (!areNonOverlapping(fragments)) {
            throw new FragmentsOverlapException();
        }
        IndigoObject restrictedContainer = restrictContainer(container, fragments);

        int[] stats = processContainer(restrictedContainer);
        int maxOpenings = stats[MAX_OPENINGS_INDEX];
        int nH          = stats[NH_INDEX];

        logger.trace("Initial molecule: " + bondsToString(restrictedContainer));
        generateMol(restrictedContainer, false, nH, maxOpenings, res, new HashSet<>());
        return res;
    }

    /**
     * Return False iff there exists two fragments in @see{fragments} with two identical elements
     * @param fragments List of fragments
     * @return False iff there exists two fragments in @see{fragments} with two identical elements
     */
    public boolean areNonOverlapping(List<IndigoObject> fragments) {
        if (fragments == null) {
            return true;
        }
        Map <IndigoObject, Set<String>> els = new HashMap<>();
        for (IndigoObject fragment : fragments) {
            els.put(fragment, new HashSet<>());
            for (IndigoObject atom : fragment.iterateAtoms()) {
                String symbol = atom.symbol();
                els.get(fragment).add(symbol);
            }
        }

        int disjunctCount = 0;

        //Test for disjunction: size of union of elements belonging to each element should be same size as sum of sizes of each fragment element sets.
        Set<String> allEls = new HashSet<>();

        for (IndigoObject frag : els.keySet()) {
            disjunctCount += els.get(frag).size();
            allEls.addAll(els.get(frag));
        }

        return disjunctCount == allEls.size();
    }

    private int[] processContainer(IndigoObject acontainer) {
        String symbol;
        int maxOpenings = 0;
        List<IndigoObject> listcont = new ArrayList<>();
        int nH = 0;
        for(IndigoObject atom: acontainer.iterateAtoms()){
            symbol = atom.symbol();
            if(symbol.equals("H")){
                nH++;
                maxOpenings--;
                listcont.add(atom);
            }
            else{
                maxOpenings += VALENCES.get(symbol);
            }
        }
        for(IndigoObject atom: listcont){
            acontainer.removeAtoms(new int[]{atom.index()});
        }
        int[] res = new int[2];
        res[MAX_OPENINGS_INDEX] = maxOpenings;
        res[NH_INDEX] = nH;
        return res;
    }

    private <K> void deltaInMap(Map<K, Integer> m, K key, int delta) {
        Integer count = m.get(key);
        if (count == null) {
            m.put(key, (count = 0));
        }
        m.put(key, count + delta);
    }


    private <K> void incrementInMap(Map<K, Integer> m, K key) {
        deltaInMap(m, key, 1);
    }

    private <K> void decrementInMap(Map<K, Integer> m, K key) {
        deltaInMap(m, key, -1);
    }

    private IndigoObject restrictContainer(IndigoObject container, List<IndigoObject> fragments) {
        IndigoObject res = indigo.createMolecule();
        logger.trace("Hydrogens count: " + container.countHydrogens());

        HashMap<String, Integer> atomsWithCount = new HashMap<>();
        for (IndigoObject atom : container.iterateAtoms()) {
            String symbol = atom.symbol();
            incrementInMap(atomsWithCount, symbol);
        }

        logger.info("Elemental composition:");
        logger.info(atomsWithCount.toString());



        if (fragments != null) {
            for (IndigoObject fragment : fragments) {
                Map<Integer, Integer> indexMappings = new HashMap<>();
                boolean[] copied = new boolean[fragment.countAtoms()];
                // First we copy all bonds from fragment
                for (IndigoObject bond : fragment.iterateBonds()) {
                    IndigoObject source = bond.source();
                    IndigoObject dest = bond.destination();
                    String sourceSymbol = source.symbol();
                    String destSymbol = dest.symbol();
                    int sourceIndex = source.index();
                    int destIndex   = dest.index();

                    copied[sourceIndex] = true;
                    copied[destIndex] = true;

                    Integer sIndex = indexMappings.get(sourceIndex);
                    if (sIndex == null) {
                        decrementInMap(atomsWithCount, sourceSymbol);
                        sIndex = res.addAtom(sourceSymbol).index();
                        indexMappings.put(sourceIndex, sIndex);
                    }

                    Integer dIndex = indexMappings.get(destIndex);
                    if (dIndex == null) {
                        decrementInMap(atomsWithCount, destSymbol);
                        dIndex = res.addAtom(destSymbol).index();
                        indexMappings.put(destIndex, dIndex);
                    }
                    addBond(res, sIndex, dIndex, bond.bondOrder());
                }

                // After this we copy all isolated atoms
                for (int i = 0; i < copied.length; i++) {
                    if (!copied[i]) {
                        String symb = fragment.getAtom(i).symbol();
                        res.addAtom(symb);
                        decrementInMap(atomsWithCount, symb);
                    }
                }
            }
        }

        // Finally we copy all atoms which are not copied before, but are in elemental composition
        for (String symb : atomsWithCount.keySet()) {
            int count = atomsWithCount.get(symb);
            if (!symb.equals("H")) {
                if (count > 0) {
                    addNAtoms(res, symb, count);
                }
            }
        }

        int hCount = container.countHydrogens() - container.countImplicitHydrogens();
        addNAtoms(res, "H", hCount);

        return res;
    }

    private void generateMol(IndigoObject acontainer,
                             boolean extraAtoms,
                             int nH,
                             int maxOpenings,
                             List<IndigoObject> accumulator,
                             Set<String> m)  {
		/*We check that the atoms are connected in the same molecules
		 * this is, they are not separated fragments*/
		/*We add hydrogens in order to check if the molecule is saturated.
		 * We will accept the molecule if the number of hydrogens necessary to saturate
		 * is the same as the hydrogens in the original formula*/
        IndigoObject acprotonate = acontainer.clone();
        int bondCount = 0;
        boolean isComplete = false;

        if(acprotonate.countHydrogens() == nH){
            isComplete = isSaturated(acprotonate);
        }

        if (isComplete && !extraAtoms) {
            if (acprotonate.countComponents() == 1) {
                    accumulator.add(acprotonate);
            }
            for (IndigoObject b : acprotonate.iterateBonds()) {
                bondCount += b.bondOrder();
            }
            if (maxOpenings > bondCount * 2) {
                generateMol(acontainer, true, nH, maxOpenings, accumulator, m);
            }
        } else {
            ArrayList<SD> extBondlist = extendMol(acontainer);

            IndigoObject molExtension;
            for (SD sourceDest : extBondlist) {
                molExtension = acontainer.clone();
                increaseBondOrder(molExtension, sourceDest.s, sourceDest.d);
                String canStr =  molExtension.canonicalSmiles();
                if (!m.contains(canStr)) {
                    m.add(canStr);
                    generateMol(molExtension, false, nH, maxOpenings, accumulator,  m);
                }
            }
        }
    }

    private IndigoObject getBond(IndigoObject mol, int a1index, int a2index) {
        for (IndigoObject bond : mol.iterateBonds()) {
            int sourceIndex = bond.source().index();
            int destIndex   = bond.destination().index();

            if (sourceIndex == a1index && destIndex == a2index ||
                    sourceIndex == a2index && destIndex == a1index) {
                return bond;
            }
        }
        return null;
    }

    private void addBond(IndigoObject mol, int source, int destination, int order) {
        mol.getAtom(source).addBond(mol.getAtom(destination), order);
    }

    private ArrayList<SD> extendMol(IndigoObject ac) {
        int vCount = ac.countAtoms();
        int [] bondCounts = new int [vCount];
        for (IndigoObject b : ac.iterateBonds()) {
            int bo = b.bondOrder();
            bondCounts[b.source().index()] += bo;
            bondCounts[b.destination().index()] += bo;
        }
        ArrayList<SD> bondList = new ArrayList<>();

        for (int i = 0; i < vCount; i++){
            IndigoObject s = ac.getAtom(i);
            if (VALENCES.get(s.symbol()) == bondCounts[i]) continue;
            for (int j = i+1; j < vCount; j++){
                IndigoObject d = ac.getAtom(j);
                if (VALENCES.get(d.symbol()) == bondCounts[j]) continue;
                bondList.add(new SD(i, j));
            }
        }
        return bondList;
    }

    private boolean isSaturated(IndigoObject mol) {
        int atomCount = mol.countAtoms();
        int[] vals = new int[atomCount];

        for (int atomIndex = 0; atomIndex < atomCount; atomIndex++) {
            vals[atomIndex]  = mol.getAtom(atomIndex).valence();
        }

        for (IndigoObject bond : mol.iterateBonds()) {
            IndigoObject source = bond.source();
            IndigoObject dest   = bond.destination();
            int bo = bond.bondOrder();

            vals[source.index()] -= bo;
            vals[dest.index()]   -= bo;
        }

        int vacant = 0;
        for (int val : vals) {
            vacant += val;
        }
        return vacant == mol.countHydrogens();
    }

    /**
     * Create IndigoObject containing atoms from given elemental composition.
     * @param ec Elemental composition
     * @return IndigoObject containing atoms in quantities given in elemental composition.
     */
    private IndigoObject loadElementalComposition(String ec) {
        String regex = "[0-9]+";
        String[] split = ec.split(regex);
        Pattern p = Pattern.compile(regex);
        Matcher m = p.matcher(ec);
        List<String> nums = new LinkedList<>();
        while (m.find()) {
            nums.add(m.group());
        }
        List<String> els = new LinkedList<>();
        els.addAll(Arrays.asList(split));
        if (nums.size() != els.size()) {
            throw new RuntimeException("Each element quantity should be specified");
        }
        IndigoObject res = indigo.createMolecule();

        int hydCount = 0;
        for (int i = 0; i < els.size(); i++) {
            String el = els.get(i);
            int cnt = Integer.valueOf(nums.get(i));
            if (!el.equals("H")) {
                addNAtoms(res, el, cnt);
            } else {
                hydCount = cnt;
            }
        }

        addNAtoms(res, "H", hydCount);
        return res;
    }

    private void addNAtoms(IndigoObject obj, String symbol, int count) {
        for (int i = 0; i < count; i++) {
            obj.addAtom(symbol);
        }
    }

    private void increaseBondOrder(IndigoObject mol, int source, int dest) {
        //Only one object bond is used, and recycled at every time we use it
        IndigoObject bondAdd = getBond(mol, source, dest);
        if(bondAdd == null){
            addBond(mol, source, dest, 1);
        } else {
            int ord = bondAdd.bondOrder();
            //TODO: to figure out why there is such limitation
            if (ord < 3) {
                bondAdd.setBondOrder(ord + 1);
            }
        }
    }

    public String bondsToString(IndigoObject mol) {
        StringBuilder res = new StringBuilder();
        for (IndigoObject bond : mol.iterateBonds()) {
            IndigoObject source = bond.source();
            IndigoObject dest   = bond.destination();
            res.append(source.symbol())
                    .append("(")
                    .append(source.index())
                    .append(")")
                    .append("-")
                    .append(dest.symbol())
                    .append("(")
                    .append(dest.index())
                    .append(")")
                    .append("[")
                    .append(bond.bondOrder())
                    .append("]");
        }
        return res.toString();
    }
}
