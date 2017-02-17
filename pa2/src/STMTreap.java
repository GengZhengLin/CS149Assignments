import java.util.concurrent.atomic.AtomicLong;

public class STMTreap implements IntSet {
    static class Node {
        final int key;
        final int priority;
        Node left;
        Node right;
        int versionNum;

        Node(final int key, final int priority) {
            this.key = key;
            this.priority = priority;
            this.versionNum = 0;
        }

        @Override
		public String toString() {
            return "Node[key=" + key + ", prio=" + priority +
                    ", left=" + (left == null ? "null" : String.valueOf(left.key)) +
                    ", right=" + (right == null ? "null" : String.valueOf(right.key)) + "]";
        }
    }
    
    static class NodePair{
    	final Node parent;
    	final Node child;
    	final int pVersion;
    	NodePair(final Node parent, final Node child){
    		this.parent = parent;
    		this.child = child;
    		if (parent != null)
    			this.pVersion = parent.versionNum;
    		else
    			this.pVersion = -1;
    	}
    }

    private AtomicLong randState = new AtomicLong();
    private Node root;

    @Override
    @org.deuce.Atomic
	public boolean contains(final int key) {
        Node node = root;
        while (node != null) {
            if (key == node.key) {
                return true;
            }
            node = key < node.key ? node.left : node.right;
        }
        return false;
    }

//    @Override
//    @org.deuce.Atomic
//	public void add(final int key) {
//        root = addImpl(root, key);
//    }
//
//    private Node addImpl(final Node node, final int key) {
//        if (node == null) {
//            return new Node(key, randPriority());
//        }
//        else if (key == node.key) {
//            // no insert needed
//            return node;
//        }
//        else if (key < node.key) {
//            node.left = addImpl(node.left, key);
//            if (node.left.priority > node.priority) {
//                return rotateRight(node);
//            }
//            return node;
//        }
//        else {
//            node.right = addImpl(node.right, key);
//            if (node.right.priority > node.priority) {
//                return rotateLeft(node);
//            }
//            return node;
//        }
//    }
    
    @org.deuce.Atomic
    NodePair addReadOnlyPhase(final int key, final int priority){
    	Node parent = null;
    	Node node = root;
    	while(node != null){
    		if (node.priority < priority){ // change to ReadWritePhase
    			return new NodePair(parent,node);
    		}
    		else if (key < node.key){
    			parent = node;
    			node = node.left;
    		}
    		else if (key > node.key){
    			parent = node;
    			node = node.right;
    		}
    		else{ // key == node.key, no need to add
    			return new NodePair(parent, node);
    		}
    	}
    	if (parent != null){
    		node = new Node(key, priority);
    		if (key < parent.key) parent.left = node;
    		else parent.right = node;
    	}
    	return new NodePair(parent,node); // node = null
    }
    
    @org.deuce.Atomic
    boolean addReadWritePhase(final NodePair nodePair, final int key, final int priority){
    	if (nodePair.parent == null){
    		root = addImpl(root, key, priority);
    	}
    	else {
    		if (nodePair.parent.versionNum != nodePair.pVersion)
    			return false;
    		else{
    			addImpl(nodePair.parent, key, priority);
//    			root = addImpl(root, key, priority);
    		}
    	}
    	return true;
    }
    
    @Override
	public void add(final int key) {
        final int priority = randPriority();
        boolean redo = false;
        do{
        	NodePair nodePair = addReadOnlyPhase(key, priority);
        	redo = !addReadWritePhase(nodePair, key, priority);
        }while(redo);
    }

    private Node addImpl(final Node node, final int key, final int priority) {
        if (node == null) {
            return new Node(key, priority);
        }
        else if (key == node.key) {
            // no insert needed
            return node;
        }
        else if (key < node.key) {
            node.left = addImpl(node.left, key, priority);
            if (node.left.priority > node.priority) {
                return rotateRight(node);
            }
            return node;
        }
        else {
            node.right = addImpl(node.right, key, priority);
            if (node.right.priority > node.priority) {
                return rotateLeft(node);
            }
            return node;
        }
    }

    private int randPriority() {
        // The constants in this 64-bit linear congruential random number
        // generator are from http://nuclear.llnl.gov/CNP/rng/rngman/node4.html
    	long v;
    	do{
    		v = randState.get();
    	}   	
        while (!randState.compareAndSet(v, v * 2862933555777941757L + 3037000493L));
        return (int)(randState.get() >> 30);
    }

    private Node rotateRight(final Node node) {
        //       node                  nL
        //     /      \             /      \
        //    nL       z     ==>   x       node
        //  /   \                         /   \
        // x   nLR                      nLR   z
        final Node nL = node.left;
        node.left = nL.right; node.versionNum += 1;
        nL.right = node; nL.versionNum += 1;
        return nL;
    }

    private Node rotateLeft(final Node node) {
        final Node nR = node.right;
        node.right = nR.left; node.versionNum += 1;
        nR.left = node; nR.versionNum += 1;
        return nR;
    }

    @org.deuce.Atomic
    NodePair removeReadOnlyPhase(int key){
    	Node parent = null;
    	Node node = root;
    	while (node != null){
    		if (key < node.key){
    			parent = node;
    			node = node.left;
    		}
    		else if (key > node.key){
    			parent = node;
    			node = node.right;
    		}
    		else{ // find node to remove
    			return new NodePair(parent, node);
    		}
    	}
    	return new NodePair(parent, node); // node = null, fail to find node to remove
    }
    
    @org.deuce.Atomic
    boolean removeReadWritePhase(NodePair nodePair, int key){
    	if (nodePair.child == null)
    		return true;
    	else if (nodePair.parent == null){
    		root = removeImpl(root,key);
    	}
    	else {
    		if (nodePair.parent.versionNum != nodePair.pVersion)
        		return false;
        	else{
        		removeImpl(nodePair.parent, key);
        	}
    	}
    	return true;
    }
    
    @Override
	public void remove(final int key) {
    	boolean redo = false;
    	do{
    		NodePair nodePair = removeReadOnlyPhase(key);
    		redo = !removeReadWritePhase(nodePair, key);
    	}while(redo);      
    }
    
//    @Override
//    @org.deuce.Atomic
//	public void remove(final int key) {
//	  	root = removeImpl(root, key);      
//    }
    
    private Node removeImpl(final Node node, final int key) {
        if (node == null) {
            // not present, nothing to do
            return null;
        }
        else if (key == node.key) {
            if (node.left == null) {
                // splice out this node
            	node.versionNum = -2;
                return node.right;
            }
            else if (node.right == null) {
            	node.versionNum = -2;
                return node.left;
            }
            else {
                // Two children, this is the hardest case.  We will pretend
                // that node has -infinite priority, move it down, then retry
                // the removal.
                if (node.left.priority > node.right.priority) {
                    // node.left needs to end up on top
                    final Node top = rotateRight(node);
                    top.right = removeImpl(top.right, key);
                    return top;
                } else {
                    final Node top = rotateLeft(node);
                    top.left = removeImpl(top.left, key);
                    return top;
                }
            }
        }
        else if (key < node.key) {
            node.left = removeImpl(node.left, key);
            return node;
        }
        else {
            node.right = removeImpl(node.right, key);
            return node;
        }
    }
}
