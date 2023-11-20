#include <iostream>
#include <cstring>

char* in() {
  char* str = new char;
  int symbol;
  int len = 0;
  symbol = getchar();
  while (!isspace(symbol)) {
    char* str1 = new char[len + 1];
    memcpy(str1, str, len);
    str1[len] = static_cast<char>(symbol);
    str = str1;
    ++len;
    symbol = getchar();
  }
  str[len] = '\0';
  return str;
}

void push(char* elem, char**& stack, int& size) {
  char** new_stack = new char* [size + 1];
  for (int i = 0; i < size; ++i) {
    new_stack[i] = stack[i];
  }
  new_stack[size] = elem;
  ++size;
  delete[] stack;
  stack = new_stack;
}

const char* pop(char**& stack, int& size) {
  if (size == 0) {
    return static_cast<const char*>("error");
  }
  char** new_stack = new char* [size - 1];
  --size;
  for (int i = 0; i < size; ++i) {
    new_stack[i] = stack[i];
  }
  char* ans = stack[size];
  delete[] stack;
  stack = new_stack;
  return ans;
}

const char* back(char** stack, int size) {
  if (size == 0) {
    return static_cast<const char*>("error");
  }
  return stack[size - 1];
}

void clear(char**& stack, int& size) {
  char** new_stack = new char* [0];
  delete[] stack;
  stack = new_stack;
  size = 0;
}


int main() {
  int size = 0;
  char** stack = new char* [0];
  while(true){
    char* query = in();
    if (strcmp(query, "push") == 0) {
      char* str = in();
      push(str, stack, size);
      std::cout << "ok" << '\n';
      continue;
    }
    if (strcmp(query,"pop") == 0) {
      std::cout << pop(stack, size) << '\n';
      continue;
    }
    if (strcmp(query,"back") == 0) {
      std::cout << back(stack, size) << '\n';
      continue;
    }
    if (strcmp(query,"size") == 0) {
      std::cout << size << '\n';
      continue;
    }
    if(strcmp(query,"clear") == 0) {
      clear(stack, size);
      std::cout << "ok" << '\n';
      continue;
    }
    if(strcmp(query, "exit") == 0) {
      delete[] stack;
      std::cout << "bye";
      break;
    }
  }
  return 0;
}
