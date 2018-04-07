// For Computer Assignment #2 of MAFS 5250 Computational Methods for Pricing Structured Products
// Date: Friday, 5 May 2017
// Revise: Saturday, 6 May 2017


public class FullyImplicitScheme {
	
	private final int factorA = 10;
	private final int factorP = 20;
	
	
	private double A0, P0;
	private int T;
	private int I, J, K, I0, J0;
	private double dA, dP, ds;
	private double r, rG;
	private double alpha, gamma, sigma;
	private double C, D, E1;
	private TridiagMatrix fullyImp;
	
	public static void main(String[] args) {
		
		// We have set all the parameters as given in Topic 4
		
		System.out.println("::::Fully Implicit Scheme for Pricing of Participating Life Insurance Policies::::\n");
		
		FullyImplicitScheme fis = new FullyImplicitScheme(100, 100, 20,
														  100, 100, 100,
														  0.05, 0.04,
														  0.3, 0.1, 0.15);
		
		System.out.println("==================== Work Element 1 ====================\n");
		
		System.out.println("A0 = 100, P0 = 100, T = 20;\n" +
		                   "I = 100, J = 100, K = 100;\n" +
				           "r = 4%, rG = 5%;\n" + 
		                   "alpha = 0.3, gamma = 0.1, sigma = 0.15.\n");
		
		System.out.println("Price = " + fis.price(true, false, true) + "\n");
		
		System.out.println("Now set (I,J) = (200, 200)");
		fis = new FullyImplicitScheme(100, 100, 20, 200, 200, 100, 0.05, 0.04, 0.3, 0.1, 0.15);
		System.out.println("Price = " + fis.price(true, false, true) + "\n");
		
		System.out.println("Now set (I,J) = (400, 400)");
		fis = new FullyImplicitScheme(100, 100, 20, 400, 400, 100, 0.05, 0.04, 0.3, 0.1, 0.15);
		System.out.println("Price = " + fis.price(true, false, true) + "\n");		
		
		System.out.println("==================== Work Element 2 ====================\n");
		
		System.out.println("Now set (I,J) = (100, 100), r = 6%, rG = 2%.");
		System.out.println("Now set (alpha, gamma) = (0.3, 0.15). Let T vary.\n");
		for (int T = 1; T <= 20; T++) {
			fis = new FullyImplicitScheme(100, 100, T, 100, 100, 100, 0.06, 0.02, 0.3, 0.1, 0.15);
			double with = fis.price(true, false, false);
			double without = fis.price(false, false, false);
			System.out.printf("T = %2d;  Value of Surrender Option = %f\n", T, with - without);
		}
		System.out.println("");
		
		System.out.println("Now set (alpha, gamma) = (0.4, 0.1). Let T vary.\n");
		for (int T = 1; T <= 20; T++) {
			fis = new FullyImplicitScheme(100, 100, T, 100, 100, 100, 0.06, 0.02, 0.4, 0.1, 0.15);
			double with = fis.price(true, false, false);
			double without = fis.price(false, false, false);
			System.out.printf("T = %2d;  Value of Surrender Option = %f\n", T, with - without);
		}
		System.out.println("");
		
		System.out.println("Now set (alpha, gamma) = (0.3, 0.15). Let T vary.\n");
		for (int T = 1; T <= 20; T++) {
			fis = new FullyImplicitScheme(100, 100, T, 100, 100, 100, 0.06, 0.02, 0.3, 0.15, 0.15);
			double with = fis.price(true, false, false);
			double without = fis.price(false, false, false);
			System.out.printf("T = %2d;  Value of Surrender Option = %f\n", T, with - without);
		}
		System.out.println("");
		
		System.out.println("==================== Work Element 3 ====================\n");

		System.out.println("Now set (alpha, gamma) = (0.3, 0.1)");
		System.out.println("Now set r = 4%. Let T vary.\n");
		for (int T = 1; T <= 20; T++) {
			fis = new FullyImplicitScheme(100, 100, T, 100, 100, 100, 0.04, 0.02, 0.3, 0.15, 0.15);
			System.out.printf("T = %2d; Price = %f\n", T, fis.price(true, false, false));
		}
		System.out.println("");
		
		System.out.println("Now set r = 7%. Let T vary.\n");
		for (int T = 1; T <= 20; T++) {
			fis = new FullyImplicitScheme(100, 100, T, 100, 100, 100, 0.07, 0.02, 0.3, 0.15, 0.15);
			System.out.printf("T = %2d; Price = %f\n", T, fis.price(true, false, false));
		}
		System.out.println("");
		
		System.out.println("Now set r = 10%. Let T vary.\n");
		for (int T = 1; T <= 20; T++) {
			fis = new FullyImplicitScheme(100, 100, T, 100, 100, 100, 0.10, 0.02, 0.3, 0.15, 0.15);
			System.out.printf("T = %2d; Price = %f\n", T, fis.price(true, false, false));
		}
		
		System.out.println("==================== FINISH ====================");

	}
	
	public FullyImplicitScheme(double A0, double P0, int T,
						 int I, int J, int K,
						 double r, double rG,
						 double alpha, double gamma, double sigma) {
		this.A0 = A0;
		this.P0 = P0;
		this.T = T;
		
		double A_ = factorA * A0;
		double P_ = factorP * P0;
		
		this.I = I;
		this.J = J;
		this.K = K;
		
		dA = A_ / I;
		dP = P_ / J;
		
		// This algorithm requires that A_, P_, I, J are chosen such that I0, J0 are integers.
		I0 = (int)(A0 / dA);
		J0 = (int)(P0 / dP);
		
		ds = 1.0 / K;
		
		this.r = r;
		this.rG = rG;
		
		this.alpha = alpha;
		this.gamma = gamma;
		this.sigma = sigma;
		
		double var = sigma * sigma;		
		
		Matrix EHG = new Matrix(I - 1, 3);
		
		C = r * ds;
		D = var * ds;
		
		for (int i = 1; i <= I - 1; i++) {

			double ii = i * i;
			
			double E =       0.5 * C * i - 0.5 * D * ii;
			double H =  1 +        C     +       D * ii;
			double G =      -0.5 * C * i - 0.5 * D * ii;
			
			if (i == 1) {
				E1 = E;
				EHG.set(i, 2, H);
				EHG.set(i, 3, G);
			}			
			else if (i < I - 1) {
				EHG.set(i, 1, E);
				EHG.set(i, 2, H);
				EHG.set(i, 3, G);
			}
			else {
				EHG.set(i, 1, E - G);
				EHG.set(i, 2, H + 2 * G);
			}
		}
		
		fullyImp = new TridiagMatrix(EHG);	
	}
	
	public double price(boolean isAmerican, boolean showProgress, boolean showTimeElapsed) {
		
		final long startTime = System.currentTimeMillis();
		
		Matrix vPlus = terminal();
		Matrix vMinus = new Matrix(I + 1, J - J0 + 1, 0, J0);
		
		for (int t = T; t > 0; t--) {
			if (showProgress)
				System.out.println("Arriving at Year " + t + "+");
			
			across(vPlus, vMinus, isAmerican);
			
			if (showProgress)
				System.out.println("Arriving at Year " + t + "-");
			
			backward(vMinus, vPlus);
			
		}		
		if (showProgress)
			System.out.println("Arriving at Year 0+");
		
		final long endTime = System.currentTimeMillis();
		
		if (showTimeElapsed)
			System.out.println("Time elapsed = " + (endTime - startTime) / 1000 + " seconds");
		
		return vPlus.get(I0, J0);
	}
	
	public Matrix terminal() {

		Matrix vPlus = new Matrix(I + 1, J - J0 + 1, 0, J0);
		
		for (int i = 0; i <= I; i++)
			for (int j = J0; j <= J; j++)
				vPlus.set(i, j, j * dP);

		return vPlus;
	}
	

	public void across(Matrix vPlus, Matrix vMinus, boolean isAmerican) {
				
		for (int i = 0; i <= I; i++)
			for (int j = J0; j <= J; j++) {
				double jTilde = jCross(i, j);
				int jFloor = Math.min((int)jTilde, J - 1);
				int jCeil = jFloor + 1; 
				vMinus.set(i, j, divisionPoint(vPlus.get(i, jFloor),
								 			   vPlus.get(i, jCeil),
								 			   jTilde - jFloor,
								 			   jCeil - jTilde));
				if (isAmerican)
					vMinus.set(i, j, Math.max(vMinus.get(i, j), j * dP));
			}
	}
	
	public double jCross(int i, int j) {
		double jTilde = j;
		double guaranteed = rG * j;
		double excess = alpha * (i * dA / dP - (gamma + 1) * j);
		jTilde += guaranteed > excess ? guaranteed : excess;
		
		return jTilde;
	}
	
	// Both linear interpolation and extrapolation can be generalized by this formula.
	// Let p be a point on the line connected by p1 and p2, whose oriented distances _
	// _ from p1 and to p2 are respectively d1 and d2.
	// We may find the coordinate of p given those of p1 and p2.
	
	public static double divisionPoint(double p1, double p2, double d1, double d2) {
		double lambda = d2 / (d1 + d2);
		return lambda * p1  + (1 - lambda) * p2; 
	}
	
	public void backward(Matrix vMinus, Matrix vPlus) {
		
		for (int j = J0; j <= J; j++) {
		
			Vector vec = vMinus.getCol(j);
			
			for (int k = 0; k <= K; k++) {
				
				Vector vec0 = vec.subVec(0, 0); vec0.align(0);
				Vector vec1 = vec.subVec(1, I - 1); vec1.align(1);
				Vector vecI = vec.subVec(I, I); vecI.align(I);
				
				vec0.set(0, vec0.get(0) * (1 - C));
				vec1.set(1, vec1.get(1) - E1 * vec0.get(0));
				vec1 = fullyImp.thomas(vec1);
				vecI.set(I, 2 * vec1.get(I - 1) - vec1.get(I - 2));
				
				vec = vec0.concat(vec1).concat(vecI);
			}
			
			vPlus.setCol(j, vec);
		}
	}
}