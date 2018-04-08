/*
 * For Assignment Bonus of MAFS 5260 Building Financial Applications with Java and VBA (Spring 2017)
 * Student Name: Fu Lei
 * Student ID: 2038 9529
 * Date: Wednesday, 3 May 2017
 */

/*
 * :::: Abstract ::::
 * This OptionPricer prices Asian options.
 * Minor changes to the structure provided in the Excel version.
 * In particular, the recursive buildTree() method in Tree Class is now a grow() method in TreeItem Class.
 * Though not required, this GUI is able to print the tree.
 * 
 */

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import javax.swing.*;
import javax.swing.border.TitledBorder;

class OptionPricerFrame extends JFrame{
	OptionPricerFrame() {
		
		// Custom classes to change font size
		class myLabel extends JLabel {
	    	myLabel(String s) {
	    		super(s);
	    		this.setFont(new Font("", Font.PLAIN, 18));
	    	}
	    }
		
	    class myTextField extends JTextField {
	    	myTextField(int n){
	    		super(n);
	    		this.setFont(new Font("", Font.PLAIN, 18));
	    	}
	    	
	    	myTextField(String s){
	    		super(s);
	    		this.setFont(new Font("", Font.PLAIN, 18));
	    	}	 
	    }
	    
	    class myRadioButton extends JRadioButton {
	    	myRadioButton(String s, boolean b) {
	    		super(s, b);
	    		this.setFont(new Font("", Font.PLAIN, 18));
	    	}	    	
	    }
	    
	    class myButton extends JButton {
	    	myButton(String s) {
	    		super(s);
	    		this.setFont(new Font("", Font.PLAIN, 18));
	    	}	    	
	    }
	    
	    class myTextArea extends JTextArea {
	    	myTextArea(int row, int column) {
	    		super(row, column);
	    		this.setFont(new Font("", Font.PLAIN, 18));
	    	}
	    }
	    
	    setTitle("Option Pricer");
		setSize(1000,1000);
	    
		// Layout the elements such as labels and textFields
		myTextField textSpotPrice = new myTextField("100");
		myTextField textRiskFreeRate = new myTextField("0.05");
		myTextField textVolatility = new myTextField("0.3");
		myTextField textMaturity = new myTextField("1");
		myTextField textStrike = new myTextField("100");
		myTextField textTimeSteps = new myTextField("3");
				
		myRadioButton callButton = new myRadioButton("Call", false);		
		myRadioButton putButton = new myRadioButton("Put", true);
		
		ButtonGroup groupOptionType = new ButtonGroup();
		groupOptionType.add(callButton);
		groupOptionType.add(putButton);
		
		JPanel panelOptionType = new JPanel();
		panelOptionType.add(callButton);
		panelOptionType.add(putButton);
		
		myRadioButton europeanButton = new myRadioButton("European", false);		
		myRadioButton americanButton = new myRadioButton("American", true);
		
		ButtonGroup groupExerciseType = new ButtonGroup();
		groupExerciseType.add(europeanButton);
		groupExerciseType.add(americanButton);
		
		JPanel panelExerciseType = new JPanel();
		panelExerciseType.add(europeanButton);
		panelExerciseType.add(americanButton);
		
		JPanel panelInput = new JPanel();
		panelInput.setLayout(new GridLayout(8, 2));
		TitledBorder borderInput = BorderFactory.createTitledBorder("Input");
		borderInput.setTitleFont(new Font("", Font.PLAIN, 18));
		panelInput.setBorder(borderInput);
		
		panelInput.add(new myLabel("Spot Price"));
		panelInput.add(textSpotPrice);
		panelInput.add(new myLabel("Risk Free Rate"));
		panelInput.add(textRiskFreeRate);
		panelInput.add(new myLabel("Volatility"));
		panelInput.add(textVolatility);
		panelInput.add(new myLabel("Option Maturity (in Years)"));
		panelInput.add(textMaturity);
		panelInput.add(new myLabel("Strike"));
		panelInput.add(textStrike);
		panelInput.add(new myLabel("Option Type"));
		panelInput.add(panelOptionType);
		panelInput.add(new myLabel("Exercise Type"));
		panelInput.add(panelExerciseType);
		panelInput.add(new myLabel("Time Steps"));
		panelInput.add(textTimeSteps);		

	    JButton buttonCalc = new myButton("Calculate");
	    
	    myTextField textDeltaT = new myTextField(20);
	    myTextField textU = new myTextField(20);
	    myTextField textD = new myTextField(20);
	    myTextField textPu = new myTextField(20);
	    myTextField textPd = new myTextField(20);
	    myTextField textDiscount = new myTextField(20);
	    myTextField textPrice = new myTextField(20);
	    
	    JPanel panelOutput = new JPanel();
	    panelOutput.setLayout(new GridLayout(7,2));
	    TitledBorder borderOutput = BorderFactory.createTitledBorder("Output");
	    borderOutput.setTitleFont(new Font("", Font.PLAIN, 18));
	    panelOutput.setBorder(borderOutput);	    
	    
	    panelOutput.add(new myLabel("Delta t"));
	    panelOutput.add(textDeltaT);
	    panelOutput.add(new myLabel("u"));
	    panelOutput.add(textU);
	    panelOutput.add(new myLabel("d"));
	    panelOutput.add(textD);
	    panelOutput.add(new myLabel("pu"));
	    panelOutput.add(textPu);
	    panelOutput.add(new myLabel("pd"));
	    panelOutput.add(textPd);
	    panelOutput.add(new myLabel("Discount"));
	    panelOutput.add(textDiscount);
	    panelOutput.add(new myLabel("Option Price"));
	    panelOutput.add(textPrice);
	    
	    JTextArea textPrint = new myTextArea(15,30);
	    JScrollPane scrollPrint = new JScrollPane(textPrint);
	    
	    JPanel panelNorth = new JPanel();
	    
	    panelNorth.setLayout(new GridLayout(1,3));
	    panelNorth.add(panelInput);
	    panelNorth.add(buttonCalc);
	    panelNorth.add(panelOutput);
	    
	    add(panelNorth, BorderLayout.NORTH);	    
	    add(scrollPrint, BorderLayout.CENTER);	    
	    
	    // ActionListenr to perform option pricing when the Calculate button is clicked
	    class CalculateListener implements ActionListener {
	        public void actionPerformed(ActionEvent event) {
	        	boolean isAmerican = americanButton.isSelected();
	        	boolean isCall = callButton.isSelected();
	        	double S0 = Double.parseDouble(textSpotPrice.getText());
	        	double K = Double.parseDouble(textStrike.getText());
	        	double r = Double.parseDouble(textRiskFreeRate.getText());
	        	double t = Double.parseDouble(textMaturity.getText());
	        	int N = Integer.parseInt(textTimeSteps.getText());
	        	double v = Double.parseDouble(textVolatility.getText());
	        	OptionPricer op = new OptionPricer(isAmerican, isCall, S0, K, r, t, N, v);
	        	textDeltaT.setText(Double.toString(op.t / op.N));
	        	textU.setText(Double.toString(op.u));
	        	textD.setText(Double.toString(1 / op.u));
	        	textPu.setText(Double.toString(op.pu));
	        	textPd.setText(Double.toString(op.pd));
	        	textDiscount.setText(Double.toString(op.D));
	        	
	        	// Re-direct console output to textPrint - idea from StackOverFlow        	
	        	textPrint.setText("");	        	
	        	// Create a stream to hold the output
	    	    ByteArrayOutputStream baos = new ByteArrayOutputStream();
	    	    PrintStream ps = new PrintStream(baos);
	    	    // IMPORTANT: Save the old System.out!
	    	    PrintStream old = System.out;
	    	    // Tell Java to use your special stream
	    	    System.setOut(ps);
	    	    // Print some output: goes to your special stream	    	    
	    	    textPrice.setText(Double.toString(op.price()));
	    	    // Put things back
	    	    System.out.flush();
	    	    System.setOut(old);
	        	textPrint.append(baos.toString() + "\n");
	        }
	     }	
	    
	    ActionListener listerner = new CalculateListener();
	    buttonCalc.addActionListener(listerner);
	    
	    pack();
	}
	
}

public class OptionPricer {
	
	boolean isAmerican, isCall;
	int callPutMultiplier;
	double S0, K, D, r, t, v, u, pu, pd;
	int N;
	
	public OptionPricer(boolean isAmerican, boolean isCall, double S0, double K, double r, double t, int N, double v) {
		// Set option type
		this.isAmerican = isAmerican;
		this.isCall = isCall;
		this.S0 = S0;
		this.K = K;
		this.N = N;
		
		callPutMultiplier = isCall ? 1 : -1;
		
		// Work out discount factor
		this.r = r;
		this.t = t;
		double dt = t / N;
		double R = Math.exp(r * dt);
		D = 1 / R;
		
		// Work out up and down probabilities
		this.v = v;
		u = Math.exp(v * Math.sqrt(dt));
		double d = 1 / u;
		pu = (R - d) / (u - d);
		pd = 1 - pu;
	}
	
	double price() {
		Tree myTree = new Tree(N);
		myTree.build(S0, u);
		evaluate(myTree.tiHead);
		myTree.print();
		return myTree.tiHead.data.optionValue;
	}
	
	
	// Recursive procedure for evaluation: left-child -> right-child -> this node
	void evaluate(TreeItem ti) {
		if (ti.leftChild == null) {
			double payoff = callPutMultiplier * (ti.data.pathData / (ti.level + 1) - K);
			ti.data.optionValue = payoff > 0 ? payoff : 0;
		}
 
		else {
			evaluate(ti.leftChild);
			evaluate(ti.rightChild);
			double binomialValue = D * (pd * ti.leftChild.data.optionValue + pu * ti.rightChild.data.optionValue);
			double exerciseValue = isAmerican ? callPutMultiplier * (ti.data.pathData / (ti.level + 1) - K) : 0;
			ti.data.optionValue = binomialValue >= exerciseValue ? binomialValue : exerciseValue;
		}
	}
	
	public static void main(String[] args) {
		JFrame frame = new OptionPricerFrame();
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setVisible(true);
	}
	
}

class Tree{
	TreeItem tiHead;
	int nLevel;
	
	Tree(int nLevel) {
		tiHead = new TreeItem();
		tiHead.level = 0;
		this.nLevel = nLevel;
	}
	
	void build(double S0, double u) {
		tiHead.grow(nLevel, S0, u);	
	}
	
	void print() {
		tiHead.printFamily();		
	}
}

class TreeItem {
	int level;
	TreeItemData data;
	TreeItem leftChild;
	TreeItem rightChild;
	TreeItem parent;
	
	TreeItem() {
		data = new TreeItemData();
	}
	
	// Allow space for left and right children
	void reproduce() {
		leftChild = new TreeItem();
		leftChild.parent = this;
		leftChild.level = this.level + 1;
		leftChild.data = new TreeItemData();
		rightChild = new TreeItem();
		rightChild.parent = this;
		rightChild.level = this.level + 1;
		rightChild.data = new TreeItemData();
	}
	
	// Recursive procedure for tree-building: this node -> left-child -> right-child
	void grow(int levelsToGrow, double price, double u) {
		data.sharePrice = price;
		this.data.pathData = price;
		if (this.parent != null)
			data.pathData += parent.data.pathData;
		
		if (levelsToGrow > 0) {
			reproduce();
			leftChild.grow(levelsToGrow - 1, price / u, u);
			rightChild.grow(levelsToGrow - 1, price * u, u);
		}
		
	}
	
	void print() {
		for (int i = 0; i < level; i++)
			System.out.print("\t");
		System.out.printf("L = %d | S = %.3f | P = %.3f | O = %.3f\n\n",
						  level, data.sharePrice, data.pathData, data.optionValue);
	}
	
	// Recursive procedure for tree-printing: left-child -> this node -> right-child
	void printFamily() {
		if (this.leftChild != null)
			this.leftChild.printFamily();
		print();
		if (this.rightChild != null)
			this.rightChild.printFamily();
	}
}

class TreeItemData {
	double sharePrice;
	double pathData;
	double optionValue;
}