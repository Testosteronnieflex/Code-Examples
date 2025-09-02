#include "LinkedList.h"

void LinkedList::push(int ID)
{
	/*Check if first element of list*/
	if (size == 0) {
		ptr_next = new Node;
		ptr_next->ID = ID;
		ptr_next->ptr_next = NULL;
		iterator = ptr_next;
		size++;
		return;
	}

	Node* next = ptr_next;
	Node* current = ptr_next;
	while (next != NULL) {
		current = next;
		next = next->ptr_next;
	}

	if (current == NULL)
		return;

	current->ptr_next = new Node;
	if (current->ptr_next != NULL) {
		current->ptr_next->ID = ID;
		current->ptr_next->ptr_next = NULL;
		size++;
	}
}

void LinkedList::pop()
{
	Node* next = ptr_next;
	Node* current = ptr_next;
	Node* prev = ptr_next;

	/*Check if zeroth element*/
	if (current->ptr_next == NULL) {
		delete current;
		size--;
		ptr_next = NULL;
		return;
	}

	while (next != NULL) {
		prev = current;
		current = next;
		next = next->ptr_next;
	}

	if (current != NULL) {
		delete current;
		size--;
		prev->ptr_next = NULL;
	}
}

int LinkedList::read(unsigned int i)
{
	unsigned int index = 0;
	Node* ptr;

	if (ptr_next == NULL) {
		printf("Empty list\n");
		return -1;
	}

	ptr = ptr_next;

	for (index; index < i; index++) {
		if (ptr->ptr_next == NULL) {
			printf("Index out of range\n");
			return -1;
		}
		ptr = ptr->ptr_next;
	}

	return ptr->ID;
}

int LinkedList::iterate()
{
	int value = -999;

	/*if (iterator == NULL) {
		iterator = ptr_next;
	}*/

	if (iterator != NULL) {
		value = iterator->ID;
		iterator = iterator->ptr_next;
	}

	return value;
}

void LinkedList::delete_element(unsigned int i)
{
	/*Check if last element*/
	if (i == (size - 1)) {
		pop();
		//printf("Delete last element\n");
		return;
	}

	/*Check if out of range*/
	if (i >= size || i < 0) {
		printf("Element out of range\n");
		return;
	}


	Node* next = ptr_next;
	Node* current = ptr_next;

	/*Check if first element*/
	if (i == 0) {
		//printf("Delete first element\n");
		next = next->ptr_next;
		if (iterator == ptr_next){ //Need to reset iterator in this case
			iterator = next;
		}
		delete ptr_next;
		ptr_next = next;
		size--;
		return;
	}

	
	for(unsigned int index = 0; index < i; index++) {
		current = next;
		next = next->ptr_next;
	}

	current->ptr_next = next->ptr_next;
	delete next;
	size--;
}

int LinkedList::find_index(int ID)
{
	int index = 0;

	Node* next = ptr_next;
	Node* current = ptr_next;
	while (next != NULL) {
		current = next;
		next = next->ptr_next;
		if (current->ID == ID) {
			return index;
		}
		index++;
	}

	return -999;
}

int LinkedList::delete_ID(int ID)
{
	int index = 0;

	Node* prev = ptr_next;
	Node* next = ptr_next;
	Node* current = ptr_next;
	while (next != NULL) {
		prev = current;
		current = next;
		next = current->ptr_next;
		if (current->ID == ID) {
			/*Check if last element*/
			if (index == (size - 1)) {
				pop();
				//printf("Delete last element\n");
				return index;
			}

			/*Check if out of range*/
			if (index >= size || index < 0) {
				printf("Element out of range\n");
				return -999;
			}

			/*Check if first element*/
			if (index == 0) {
				//printf("Delete first element\n");
				if (iterator == ptr_next) { //Need to reset iterator in this case
					iterator = next;
				}
				delete ptr_next;
				ptr_next = next;
				size--;
				return index;
			}

			prev->ptr_next = next;
			delete current;
			size--;
			return index;
		}
		index++;
	}

	return -999;
}


void LinkedList::empty()
{
	Node* next = ptr_next;
	Node* current = ptr_next;
	while (next != NULL) {
		current = next;
		next = current->ptr_next;
		delete current;
	}

	size = 0;
	ptr_next = NULL;
}

Cells::Cells()
{
	nx = 0;
	ny = 0;
	nz = 0;
	sx = 0.0;
	sy = 0.0;
	sz = 0.0;
	cell_lists = NULL;
}

Cells::Cells(int Nx, int Ny, int Nz)
{
	nx = Nx;
	ny = Ny;
	nz = Nz;
	sx = 0.0;
	sy = 0.0;
	sz = 0.0;

	cell_lists = new LinkedList[Nx * Ny * Nz];
}

void Cells::create(int Nx, int Ny, int Nz)
{
	nx = Nx;
	ny = Ny;
	nz = Nz;

	cell_lists = new LinkedList[nx * ny * nz];
}

LinkedList* Cells::access(int ix, int iy, int iz)
{
	if (ix > nx - 1 || iy > ny - 1 || iz > nz - 1 || ix < 0 || iy < 0 || iz < 0) {
		printf("Arr3D index out of range! \n");
		return NULL;
	}
	int index = ix * ny * nz + iy * nz + iz;

	return &cell_lists[index];
}