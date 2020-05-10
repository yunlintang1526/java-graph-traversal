/*
 * Name: Yunlin Tang
 * PID: a14664383
 */

import java.util.ArrayList;
import java.util.PriorityQueue;

/**
 * the helping class to construct the Graph class
 *
 * @author Yunlin Tang
 * @since ${12.2}
 */
public class Vertex {

    private final String name; // the name of this vertex
    private final int x; // the x coordinates of this vertex on map
    private final int y; // the y coordinates of this vertex on map

    private ArrayList<Edge> nextVertex; // list store all the vertex connect
                                        // to the current vertex
    private Vertex prev; // the previous vertex which follows the path


    /**
     * constructor of Vertex class
     * @param name the name of vertex
     * @param x the x-coordinate
     * @param y the y-coordinate
     */
    public Vertex(String name, int x, int y) {
        this.name = name;
        this.x = x;
        this.y = y;
        nextVertex = new ArrayList<>(); // initialize the children vertex list
    }

    /**
     * getter for name
     * @return the name of this vertex
     */
    public String getName() {
        return name;
    }

    /**
     * getter for x-coordinate
     * @return the x coordinate for this vertex
     */
    public int getX() {
        return x;
    }

    /**
     * the getter method of y-variable
     * @return the y-coordinate for this vertex
     */
    public int getY() {
        return y;
    }

    /**
     * setter for next arrayList
     * @param v the children vertex
     * @param weight weight on the edge
     */
    public void setNext(Vertex v, Double weight){
        Edge e = new Edge(this, v, weight); // construct a new Edge
        nextVertex.add(e); // add the edge to the list
    }

    /**
     * getter for the next list
     * @return return the next list which includes the children vertex
     */
    public ArrayList<Edge> getNext(){
        return nextVertex;
    }

    /**
     * setter method for previous vertex
     * @param v the input vertex to be set as the previous
     */
    public void setPrev(Vertex v){
        prev = v;
    }

    /**
     * the getter method for previous vertex
     * @return the previous vertex
     */
    public Vertex getPrev(){
        return prev;
    }

    /**
     * create a hashcode for this vertex
     * @return the hash code of vertex
     */
    @Override
    public int hashCode() {
        // we assume that each vertex has a unique name
        return name.hashCode();
    }

    /**
     * method to determine if this vertex equals to other vertex
     * @param o other vertex
     * @return the result if they are equal or not
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null) {
            return false;
        }
        if (!(o instanceof Vertex)) {
            return false;
        }
        Vertex oVertex = (Vertex) o;

        return name.equals(oVertex.name) && x == oVertex.x && y == oVertex.y;
    }

    /**
     * convert the vertex to string representation
     * @return the string representation
     */
    public String toString() {
        return name + " (" + x + ", " + y + ")";
    }

}