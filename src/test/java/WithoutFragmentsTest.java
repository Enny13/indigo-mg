import com.ggasoftware.indigo.Indigo;
import org.apache.log4j.Logger;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import org.indigomg.FragmentsOverlapException;
import org.indigomg.MoleculesGenerator;

import java.util.HashMap;
import java.util.Map;

/**
 * Tests generating molecules without 'fragments' restriction
 * Created by Artem Malykh on 18.09.15.
 */
public class WithoutFragmentsTest {
    private static Map<String, Long> elCompositionsWithExpectedResults = new HashMap<>();
    private static Indigo indigo;
    private static Logger logger = Logger.getLogger(WithoutFragmentsTest.class);

    static {
        elCompositionsWithExpectedResults.put("C3H6",    2L);
        elCompositionsWithExpectedResults.put("C6H8",    159L);
        elCompositionsWithExpectedResults.put("C3H4O3",  152L);
        elCompositionsWithExpectedResults.put("C2H5N1O2", 84L);
    }

    @BeforeClass
    public static void before() {
        indigo = new Indigo();
    }

    @Test
    public void test() throws FragmentsOverlapException {
        MoleculesGenerator gen = new MoleculesGenerator(indigo);

        for (Map.Entry<String, Long> compWithRes : elCompositionsWithExpectedResults.entrySet()) {
            String  comp = compWithRes.getKey();
            long res  = compWithRes.getValue();
            logger.info("Composition: " + comp);
            long before = System.currentTimeMillis();
            long count = gen.generateMolecules(comp, (String)null).size();
            logger.info("Total molecules: " + count + ". Took " + (System.currentTimeMillis() - before) + " millis.");
            Assert.assertEquals(res, count);
        }
    }
}
