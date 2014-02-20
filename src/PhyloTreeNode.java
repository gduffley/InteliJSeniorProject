/**
 * Created by Gordon on 2/11/14.
 */
import java.util.ArrayList;

public class PhyloTreeNode {
    private String sequence;
    private String name;
    private PhyloTreeNode parent;
    private ArrayList<PhyloTreeNode> children = new ArrayList<PhyloTreeNode>();

    public PhyloTreeNode(String name, String sequence){
        this.name = name;
        this.sequence = sequence;
        this.parent = null;
    }

    public void setName(String name){
        this.name = name;
    }

    public void setParent(PhyloTreeNode parent){
        this.parent = parent;
    }

    public void addChild(PhyloTreeNode a){
        children.add(a);
        a.setParent(this);
    }

    public void removeChild(PhyloTreeNode a){
        this.children.remove(a);
    }

    public void setSequence(String sequence){
        this.sequence = sequence;
    }

    public String getName(){
        return this.name;
    }

    public String getSequence(){
        return this.sequence;
    }

    public PhyloTreeNode getParent() {return this.parent;}

    public ArrayList<PhyloTreeNode> getChildren() {return this.children;}
}

