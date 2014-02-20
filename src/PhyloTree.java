/**
 * Created by Gordon on 2/11/14.
 */

public class PhyloTree {
    private PhyloTreeNode root;


    public PhyloTree(String name, String sequence){
        root = new PhyloTreeNode(name, sequence);
    }

    public void addNode(PhyloTreeNode newRoot){
        newRoot.addChild(this.root);
        root = newRoot;
    }
    public PhyloTreeNode getRoot(){
        return this.root;
    }


}