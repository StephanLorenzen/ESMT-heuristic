#ifndef UTILS_HEAP_H
#define UTILS_HEAP_H

#include <vector>

/*
 * @namespace Utils
 */
namespace Utils {
  /**
   * @class Heap
   * This class implements a heap data structure, as described
   * in T. H. Cormens et al. (Introduction to Algorithms).
   */
  template <class T>
  class Heap {
  public:
    /**
     * Constructs an empty heap.
     *
     * @param hFunc  The comparision function, defining the heap
     */
    Heap(bool (* hFunc)(const T&, const T&));
    /**
     * Constructs the heap from the given (unordered) list.
     *
     * @param hFunc  The comparision function, defining the heap
     * @param input  The initial elements of the heap, in an unordered list
     */
    Heap(bool (* hFunc)(const T&, const T&), std::vector<T> &input);
    /**
     * Destructor
     */
    ~Heap();

    /**
     * Inserts a value into the heap.
     * @param val  The value to be inserted.
     */
    void insert(T &val);
    
    /**
     * Returns the root element of the heap, without removing it.
     *
     * @return The root element of the heap.
     */
    T &peek();

    /**
     * Removes and returns the root element of the heap
     *
     * @return The root element of the heap.
     */
    T &extract();
    
    /**
     * Gets the size (number of elements) of the heap.
     *
     * @return The size of the heap
     */
    unsigned int size() const;
  protected:
  private:
    
    /**
     * Helper function used for reordering the heap in order to satisfy the heap invariant. 
     */
    void heapify(unsigned int i);
    
    /**
     * Getter for the parent of i
     */
    unsigned int parent(unsigned int i);
    /**
     * Getter for the left child of i
     */
    unsigned int left(unsigned int i);
    /**
     * Getter for the right child of i
     */
    unsigned int right(unsigned int i);
    
    /** The comparison function */
    bool (* hFunc)(const T&, const T&);
    /** Pointers to the values stored in this heap. */
    std::vector<T*> vals;
  };

}

#endif // UTILS_HEAP_H
