public class Matrix {
	
	// m = number of rows
	// n = number of columns
	// i = index of a row
	// j = index of a column
	
	int m, n;
	
	// We allow row and column indices to start from, say, 0.
	// By default, they start from one.
	
	// Auto boxing and unboxing may occur between int and Integer
	
	Integer firstRowIndex;
	Integer firstColIndex;
	Integer lastRowIndex;
	Integer lastColIndex;
	
	double[][] A;
	
	public Matrix(int m, int n) {
		
		this(m, n, 1, 1);
	}
	
	public Matrix(int m, int n, int firstRowIndex, int firstColIndex) {
		
		this.m = m;
		this.n = n;
		A = new double[m][n];
		
		align(firstRowIndex, firstColIndex);
	}
	
	public void align(int firstRowIndex, int firstColIndex) {

		this.firstRowIndex = firstRowIndex;
		this.firstColIndex = firstColIndex;
		
		lastRowIndex = firstRowIndex + m - 1;
		lastColIndex = firstColIndex + n - 1;		
	}
	
	public void align(Matrix mat) {
		
		align(mat.firstRowIndex, mat.firstColIndex);
	}
	
	public double get(int i, int j) {
		
		return A[i - firstRowIndex][j - firstColIndex];
	}
	
	public Vector getRow(int i) {
		
		Vector vec = new Vector(n, firstRowIndex);
		
		for (int j = firstColIndex; j <= lastColIndex; j++)
			vec.set(j, get(i,j));
		
		return vec;
	}
	
	public Vector getCol(int j) {
		
		Vector vec = new Vector(m, firstRowIndex);
		for (int i = firstRowIndex; i <=  lastRowIndex; i++)
			vec.set(i, get(i,j));
		
		return vec;
	}
	
	public void set(int i, int j, double x) {
	
		A[i - firstRowIndex][j - firstColIndex] = x;
	}
	
	public void setCol(int j, Vector vec) {
		
		int vecFirstIndex = vec.firstIndex;
		
		vec.align(this);
		
		for (int i = firstRowIndex; i <= lastRowIndex; i++)
			set(i, j, vec.get(i));
		
		vec.align(vecFirstIndex);
	}
	
	public void rowMul(int i, double c) {
		
		for (int j = firstColIndex; j <= lastColIndex; j++)
			set(i, j, c * get(i,j));
	}
	
	public void rowDiv(int i, double c) {
		
		rowMul(i, 1 / c);
	}
	
	public void rowAdd(int i, int a, double c) {
		
		// Add (c * Row a) to Row i
		
		for (int j = firstColIndex; j <= lastColIndex; j++)
			set(i, j, get(i,j) + c * get(a,j)); 
	}
	
	public void rowSub(int i, int a, double c) {
		
		rowAdd(i, a, -1 * c);
	}
	
	public Matrix mul(Matrix mat) {
		
		if (mat.m != n)
			return null;
		
		int matFirstRowIndex = mat.firstRowIndex;
		int matFirstColIndex = mat.firstColIndex;
		
		mat.align(this);
		
		Matrix prod = new Matrix(m, mat.n, firstRowIndex, firstColIndex);
		for (int i = firstRowIndex; i <= lastRowIndex; i++)
			for(int j = firstColIndex; j <= mat.lastColIndex; j++)
				prod.set(i, j, getRow(i).dotMul(mat.getCol(j)));
		
		mat.align(matFirstRowIndex, matFirstColIndex);
		
		return prod;
	}
	
	public Matrix augment(Vector beta) {
		
		// Return an augmented matrix [ A | beta ]
		
		if (beta.m != m)
			return null;
		
		int vecFirstIndex = beta.firstIndex;
		
		beta.align(this);
		
		Matrix aug  =  new Matrix(m, n + 1, firstRowIndex, firstColIndex);
		for (int i = aug.firstRowIndex; i <= aug.lastRowIndex; i++)
			for (int j = aug.firstColIndex; j <= aug.lastColIndex; j++)
				aug.set(i, j, j <= lastColIndex ? get(i,j) : beta.get(i));
		
		beta.align(vecFirstIndex);
		
		return aug;
	}
	
	public void print() {
		print(true);
	}
	
	public void print(boolean showIndices) {
		
		String tab = "    ";
		String newline = "\n";
		String nothing = "";
		
		if (showIndices) {
			System.out.print(tab + tab);
			for (int j = firstColIndex; j <= lastColIndex; j++)
				System.out.printf("[%2d]%s", j, j == lastColIndex ? newline : tab);
		}
		
		
		for (int i = firstRowIndex; i <= lastRowIndex; i++) 
			for (int j = firstColIndex; j <= lastColIndex;  j++) {
				if (showIndices && j == firstColIndex)
					System.out.printf("[%2d]", i);
				System.out.printf("%8.1f%s", get(i,j), j == lastColIndex ? newline : nothing);
			}
	}
}

class SquareMatrix extends Matrix {

	public SquareMatrix(int n) {
		super(n, n);
	}
	
	public SquareMatrix(int n, int firstRowIndex, int firstColIndex) {
		super(n, n, firstRowIndex, firstColIndex);
	}
}

class Triple {
	
	double a, b, c;
	
	Triple(double a, double b, double c) {
		this.a = a;
		this.b = b;
		this.c = c;
	}	
}

class TridiagMatrix extends SquareMatrix {

	public TridiagMatrix(int n) {
		
		super(n);
	}
	
	public TridiagMatrix(int n, Triple t) {
		
		super(n);
		
		for (int i = firstRowIndex; i <= lastRowIndex; i++) {
			if (i > firstRowIndex)
				set(i, i - 1, t.a);
			set(i, i, t.b);
			if (i < lastRowIndex)
				set(i, i + 1, t.c);
		}
	}
	
	public TridiagMatrix(Matrix mat) {
		
		// Construct a tridiagonal matrix from a m-by-3 matrix
				
		super(mat.m, mat.firstRowIndex, mat.firstColIndex);
		
		for (int i = firstRowIndex; i <= lastRowIndex; i++) {
			if (i > firstRowIndex)
				set(i, i - 1, mat.get(i, mat.firstColIndex));
			set(i, i, mat.get(i, mat.firstColIndex + 1));
			if (i < lastRowIndex)
				set(i, i + 1, mat.get(i, mat.firstColIndex + 2));
		}
		
	}
	
	public Vector thomas(Vector beta) {
		
		if (beta.m != m)
			return null;

		Matrix aug = augment(beta);
		
		// Solve Ax = beta for x
		aug.rowDiv(firstRowIndex, aug.get(firstRowIndex, firstColIndex));
		
		// Forward
		for (int i = firstRowIndex + 1; i <= lastRowIndex; i++) {
			aug.rowSub(i, i - 1, aug.get(i, i - 1));
			aug.rowDiv(i, aug.get(i,i));
		}
		
		// Backward
		for (int i = lastRowIndex - 1; i >= firstRowIndex; i--) {
			aug.rowSub(i, i + 1, aug.get(i, i + 1));
		}
		
		return aug.getCol(aug.lastColIndex);
	}
}

class Vector extends Matrix {

	// We see Vector as a m-by-1 Matrix.
	// We allow its firstRowIndex to vary, but require its firstColIndex to be 1.
	
	Integer firstIndex = firstRowIndex;
	Integer lastIndex = lastRowIndex;
	
	public Vector(int m) {
		
		this(m, 1);
	}
	
	public Vector(int m, int firstIndex) {
		
		super(m, 1, firstIndex, 1);
	}
	
	public void align(int firstIndex) {
		
		align(firstIndex, 1);
	}
	
	public void align(Vector vec) {
		
		align(vec.firstIndex);
	}
	
	public void align(Matrix mat) {
		
		align(mat.firstRowIndex);
	}
	
	public int length() {
		
		return m;
	}
	
	public double get(int i) {
		
		return get(i, 1);
	}
	
	public void set(int i, double x) {
		
		set(i, 1, x);
	}
	
	public void set(Vector vec) {
		
		int vecFirstIndex = firstIndex;

		vec.align(this);
		
		for (int i = firstRowIndex; i <= lastRowIndex; i++)
			set(i, vec.get(i));
		
		vec.align(vecFirstIndex);		
	}
	
	public Vector subVec(int a, int b) {
		
		// Generate a sub-vector from a to b
		Vector sub = new Vector(b - a + 1, a);
		
		for (int i = a; i <= b; i++)
			sub.set(i, get(i));
		
		sub.align(this);
		
		return sub;
	}
	
	public Vector concat(Vector vec) {
		
		int totLen = m + vec.length();
		
		Vector tot = new Vector(totLen, firstIndex);
		
		int vecFirstIndex = vec.firstIndex;
		
		vec.align(this);
		
		for (int i = tot.firstIndex; i <= lastIndex; i++)
			tot.set(i, get(i));
		for (int i = lastIndex + 1; i <= tot.lastIndex; i++)
			tot.set(i, vec.get(i - m));
		
		vec.align(vecFirstIndex);
		
		return tot;
	}
	
	
	public double dotMul(Vector vec) {
		
		int vecFirstIndex = vec.firstIndex;
		
		vec.align(this);
		
		double sum = 0;
		for (int i = firstIndex; i <= lastIndex; i++)
			sum += get(i) * vec.get(i);
		
		vec.align(vecFirstIndex);
		
		return sum;
	}
}