#include <cstdlib>
#include <vector>
#include <math.h>

#include "steiner/steiner_tree.hpp"
#include "steiner/utils/heap.hpp"

/*
 * Implementation of Heap::Heap(...) for constructing empty heap
 */
template <class T>
Utils::Heap<T>::Heap(bool (* hFunc)(const T&, const T&))
  : hFunc(hFunc)
{ }

/*
 * Implementation of Heap::Heap(...) for constructing heap with initial values
 */
template <class T>
Utils::Heap<T>::Heap(bool (* hFunc)(const T&, const T&), std::vector<T> &input)
  : hFunc(hFunc), vals(input.size(), NULL)
{
  for(unsigned int i=0; i<input.size(); i++)
    this->vals[i] = &input[i];
  for(int i = floor(input.size() / 2)-1; i >= 0; i--) {
    this->heapify(i);
  }
}

/*
 * Implementation of Heap::~Heap()
 */
template <class T>
Utils::Heap<T>::~Heap() { }

/*
 * Implementation of Heap::insert(...)
 */
template <class T>
void Utils::Heap<T>::insert(T &val) {
  this->vals.push_back(&val);
  unsigned int i = this->size()-1;
  while(i > 0 && this->hFunc(*this->vals[i],*this->vals[this->parent(i)])) {
    T* tmp = this->vals[i];
    this->vals[i] = this->vals[this->parent(i)];
    this->vals[this->parent(i)] = tmp;
    i = this->parent(i);
  }
}

/*
 * Implementation of Heap::peek()
 */
template <class T>
T &Utils::Heap<T>::peek() {
  if(this->size() < 1)
    throw "Heap underflow";
  return *this->vals[0];
}

/*
 * Implementation of Heap::extract()
 */
template <class T>
T &Utils::Heap<T>::extract() {
  if(this->size() < 1)
    throw "Heap underflow";
  T &val = this->peek();
  this->vals[0] = this->vals.back();
  this->vals.pop_back();
  if(this->size() > 1)
    this->heapify(0);
  return val;
}

/*
 * Implementation of Heap::size()
 */
template <class T>
unsigned int Utils::Heap<T>::size() const {
  return this->vals.size();
}

/*
 * Implementation of Heap::heapify(...)
 */
template <class T>
void Utils::Heap<T>::heapify(unsigned int i) {
  unsigned int l, r, m;
  l = this->left(i);
  r = this->right(i);
  m = l < this->size() && this->hFunc(*this->vals[l], *this->vals[i]) ? l : i;
  if(r < this->size() && this->hFunc(*this->vals[r], *this->vals[m]))
    m = r;
  if(m != i) {
    // Exchange m with n
    T* tmp = this->vals[m];
    this->vals[m] = this->vals[i];
    this->vals[i] = tmp;
    this->heapify(m);
  }
}

/*
 * Implementation of Heap::parent(...)
 */
template <class T>
unsigned int Utils::Heap<T>::parent(unsigned int i) {
  return floor((i-1)/2);
}

/*
 * Implementation of Heap::left(...)
 */
template <class T>
unsigned int Utils::Heap<T>::left(unsigned int i) {
  return 2*i+1;
}

/*
 * Implementation of Heap::right(...)
 */
template <class T>
unsigned int Utils::Heap<T>::right(unsigned int i) {
  return 2*i+2;
}

template class Utils::Heap<SteinerTree>;
