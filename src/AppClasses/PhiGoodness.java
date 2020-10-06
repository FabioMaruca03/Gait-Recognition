package AppClasses;

public class PhiGoodness
{
	private double phi;
	private boolean rightClassified;
	
	public PhiGoodness(double phi, boolean rightClassified)
	{
		this.setPhi(phi);
		this.setRightClassified(rightClassified);
	}

	public double getPhi() {return phi;}
	public void setPhi(double phi) {this.phi = phi;}
	public boolean isRightClassified() {return rightClassified;}
	public void setRightClassified(boolean rightClassified) {this.rightClassified = rightClassified;}
}