package AppClasses;

public class MaximumPos implements Comparable<MaximumPos>
{
	private int value;
	private int pos;
	
	public MaximumPos(int v, int p)
	{
		setValue(v);
		pos=p;
	}

	@Override
	public int compareTo(MaximumPos arg0)
	{
		if(getValue()>arg0.getValue())
			return 1;
		return (getValue()==arg0.getValue()) ? 0:-1;
	}

	public int getValue() {return value;}
	public void setValue(int value) {this.value = value;}

	public int getPos() {return pos;}
}
