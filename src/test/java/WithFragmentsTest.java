import com.ggasoftware.indigo.Indigo;
import com.ggasoftware.indigo.IndigoObject;
import org.apache.log4j.Logger;
import org.junit.Assert;
import org.junit.Test;
import org.indigomg.FragmentsOverlapException;
import org.indigomg.MoleculesGenerator;

import java.util.*;

/**
 * Tests generation molecules with fragments restriction.
 * Created by Artem Malykh on 18.09.15.
 */
public class WithFragmentsTest {
    private static Map<CompWithRestrictions, Long> elCompositionsWithExpectedResults = new HashMap<>();
    private static IndigoObject CO2;
    private static Indigo indigo = new Indigo();
    private static Logger logger = Logger.getLogger(WithoutFragmentsTest.class);

    static {
        CO2 = indigo.createMolecule();
        IndigoObject c = CO2.addAtom("C");
        IndigoObject o1 = CO2.addAtom("O");
        IndigoObject o2 = CO2.addAtom("O");

        c.addBond(o2, 2);
        c.addBond(o1, 1);
    }

    private static class CompWithRestrictions {
        private String composition;
        private List<IndigoObject> restrictions;

        public CompWithRestrictions(String composition, List<IndigoObject> restrictions) {
            this.composition = composition;
            this.restrictions = restrictions;
        }
    }

    static {
        elCompositionsWithExpectedResults.put(new CompWithRestrictions("N1H5C2O2", Arrays.asList(new IndigoObject[] {CO2})), 6L);
    }

    @Test
    public void test() throws FragmentsOverlapException {
        MoleculesGenerator gen = new MoleculesGenerator(indigo);

        for (Map.Entry<CompWithRestrictions, Long> compWithRes : elCompositionsWithExpectedResults.entrySet()) {
            CompWithRestrictions  comp = compWithRes.getKey();
            long expectedCount  = compWithRes.getValue();
            logger.info("Composition: " + comp);
            logger.info("Should contain following substructures: ");
            for (IndigoObject restriction : comp.restrictions) {
                logger.info(restriction.smiles());
            }

            long before = System.currentTimeMillis();
            long count = gen.generateMolecules(comp.composition, comp.restrictions).size();
            logger.info("Total molecules: " + count + ". Took " + (System.currentTimeMillis() - before) + " millis.");
            Assert.assertEquals(expectedCount, count);
        }
    }
}
