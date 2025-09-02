#pragma once
#include <stdio.h>

struct Node{
	int ID;
	Node* ptr_next;
};

class LinkedList{
	unsigned int size;
	Node* ptr_next;
	Node* iterator;

public:
	/*--------------Constructor-----------------*/
	LinkedList() { size = 0; ptr_next = NULL; iterator = NULL; }


	/*Returns length of LinkedList*/
	unsigned int len() { return size; }

	/*Appends node with value ID to end of list*/
	void push(int ID);

	/*Remove last item of list*/
	void pop();

	/*Read value from index i*/
	int read(unsigned int i);

	/*Gives next value in the linked list eacht time this is called until end of list
	note: returns -999 if a value could not be read (e.g. iterated beyond end of list)*/
	int iterate();

	void reset_iterate() { iterator = ptr_next; }

	/*Delete element with index i*/
	void delete_element(unsigned int i);

	/*Finds index of element containing the ID, returns -999 if element not in list*/
	int find_index(int ID);

	/*Empty the list of all its elements*/
	void empty();
};

