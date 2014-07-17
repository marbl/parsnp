/////////////////////////////////////////
// Converter.h
// simple Stack for MUM calculator
/////////////////////////////////////////

// See the LICENSE file included with this software for license information.

#ifndef STACK_H
#define STACK_H


#include <stdlib.h>
#include <assert.h>
#include <iostream>

//the following 2 lines are necessary for
//templating the stack
template <class T>


class Stack
{

public:
	Stack(void);
	~Stack(void);
	void push(T);
	T pop(void);
	//this is const since we dont change 
	//anything in the stack w/ peek
	T peek(void) const;
	bool full(void)const;//test if stack is full
	bool empty(void)const;//test is stack is empty

private:

	T *arr; //templated pointer to stack
	int top; //holds position of element to be popped
	int size; //size of stack


};

template <class T>
Stack<T>::Stack(void)
{
  size = 80;
  arr = new T[size];
  assert(arr);
  top = 0;

}

//bool Stack<T>::full(void)const
//determines if stack is empty
//input:  top
//output:  nothing
//returns:  true if empty, false otherwise
//
template <class T>
bool Stack<T>::empty(void)const
{
  return(top <= 0);
}

//bool Stack<T>::full(void)const
//determines if stack is full
//input: top, size
//output: nothing
//returns: true if full, false otherwise
//
template <class T>
bool Stack<T>::full(void)const
{
  return(top >= size);
}


//void Stack<T>::push(T toma)
//pushes an element of type T onto the stack
//increases the size of the stack if full
//increase the top by 1
//input: arr, size, toma, top
//output: the new values of arr, size, top
//returns: nothing
//
template <class T>
void Stack<T>::push(T toma)
{
  if(full())
  {
    int nsize = size + 5;
    T *narr = new T[nsize];
    assert(narr);
    for(int i = 0; i < size; i++)
    {
      narr[i] = arr[i];
    }
    delete [] arr;
    arr = narr;
    size = nsize;
  }
  arr[top] = toma;
  top++;
}

//T Stack<T>::pop(void)
//pops an element of type T from the stack
//gives an error message if stack is empty
//decreases the top by 1
//input: arr, top
//output: arr & top modified
//returns: the element of type T from the Stack
//
template <class T>
T Stack<T>::pop(void)
{
  if(empty())
  {
    //cerr << "Can't pop an empty stack!" << endl;
    exit(1);
  }
  return arr[--top];
}

//T Stack<T>::peek(void)const
//peeks at an element of type T from the stack
//gives an error message if stack is empty
//does not decrease the top by 1, peek simply
//looks at the element, doesnt modify the stack
//in any way
//input:  arr, top
//output: nothing
//returns;  the element of type T from the Stack
//
template <class T>
T Stack<T>::peek(void)const
{
  if(empty())
  {
    //cerr << "Can't peek at an empty stack!" << endl;
    exit(1);
  }

  //peek doesnt delete the element looked at
  return arr[top-1];
}

//Stack<T>::~Stack(void)
//destructor for Class Stack
//deletes the memory for the stack to
//avoid a memory leak
//input:  arr
//output: arr deleted
//returns: nothing
//
template <class T>
Stack<T>::~Stack(void)
{
  delete [] arr;
}

#endif
  
