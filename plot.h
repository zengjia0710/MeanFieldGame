#pragma once

#include<Python.h>
#include<string>
using namespace std;
#include <iostream>

/*to plot through Python*/
template<class T>
string arr_to_string_list(T* arr, int N) {
	string s = "[";
	for (int i = 0; i < N; ++i) {
		s += to_string(arr[i]);
		if (i != N - 1) s += ",";
	}
	s += "]";
	return s;
}

/*one line*/
template<class T, class V = int>
void plot(T* x, int N1, string title, V* y = NULL, bool equal = false) {
	PyRun_SimpleString("import matplotlib.pyplot as plt");
	if (equal) {
		PyRun_SimpleString("plt.axis(\"equal\")");
	}

	string cmd = "plt.plot(";
	string s1 = arr_to_string_list(x, N1);
	if (y != NULL) {
		string s2 = arr_to_string_list(y, N1);
		cmd += (s1 + "," + s2 + ")");
		PyRun_SimpleString(cmd.c_str());
	}
	else {
		cmd += (s1 + ")");
		PyRun_SimpleString(cmd.c_str());
	}
	string cmd2 = "plt.title('";
	cmd2 += (title + "')");
	PyRun_SimpleString(cmd2.c_str());
	PyRun_SimpleString("plt.xlabel('time t')");
	PyRun_SimpleString("plt.show()");
}

/*two lines*/
template<class T, class V = int>
void plot(T* x, int N1, V* y, string label1, string title, V* z, string label2, bool equal = false) {
	PyRun_SimpleString("import matplotlib.pyplot as plt");
	if (equal) {
		PyRun_SimpleString("plt.axis(\"equal\")");
	}

	string cmd = "plt.plot(";
	string s1 = arr_to_string_list(x, N1);

	string s2 = arr_to_string_list(y, N1);
	cmd += (s1 + "," + s2 + ", label='" + label1 + "')");
	PyRun_SimpleString(cmd.c_str());

	if (z != NULL) {
		string cmd1 = "plt.plot(";
		string s3 = arr_to_string_list(x, N1);

		string s4 = arr_to_string_list(z, N1);
		cmd1 += (s3 + "," + s4 + ", label='" + label2 + "')");
		PyRun_SimpleString(cmd1.c_str());
	}

	string cmd2 = "plt.title('";
	cmd2 += (title + "')");
	PyRun_SimpleString(cmd2.c_str());
	PyRun_SimpleString("plt.xlabel('time t')");
	PyRun_SimpleString("plt.legend()");


	PyRun_SimpleString("plt.show()");
}

/*three lines*/
template<class T, class V = int>
void plot(T* x, int N1, V* y, string label1, string title, V* z, string label2, V* u, string label3, bool equal = false) {
	PyRun_SimpleString("import matplotlib.pyplot as plt");
	if (equal) {
		PyRun_SimpleString("plt.axis(\"equal\")");
	}

	string cmd = "plt.plot(";
	string s1 = arr_to_string_list(x, N1);

	string s2 = arr_to_string_list(y, N1);
	cmd += (s1 + "," + s2 + ", label='" + label1 + "')");
	PyRun_SimpleString(cmd.c_str());

	if (z != NULL) {
		string cmd1 = "plt.plot(";
		string s3 = arr_to_string_list(x, N1);

		string s4 = arr_to_string_list(z, N1);
		cmd1 += (s3 + "," + s4 + ", label='" + label2 + "')");
		PyRun_SimpleString(cmd1.c_str());
	}

	if (u != NULL) {
		string cmd2 = "plt.plot(";
		string s5 = arr_to_string_list(x, N1);

		string s6 = arr_to_string_list(u, N1);
		cmd2 += (s5 + "," + s6 + ", label='" + label3 + "')");
		PyRun_SimpleString(cmd2.c_str());
	}

	string cmd3 = "plt.title('";
	cmd3 += (title + "')");
	PyRun_SimpleString(cmd3.c_str());
	PyRun_SimpleString("plt.xlabel('time t')");
	PyRun_SimpleString("plt.legend()");


	PyRun_SimpleString("plt.show()");
}

/*four lines*/
template<class T, class V = int>
void plot(T* x, int N1, V* y, string label1, string title, V* z, string label2, V* u, string label3, V* v, string label4, bool equal = false) {
	PyRun_SimpleString("import matplotlib.pyplot as plt");
	if (equal) {
		PyRun_SimpleString("plt.axis(\"equal\")");
	}

	string cmd = "plt.plot(";
	string s1 = arr_to_string_list(x, N1);

	string s2 = arr_to_string_list(y, N1);
	cmd += (s1 + "," + s2 + ", label='" + label1 + "')");
	PyRun_SimpleString(cmd.c_str());

	if (z != NULL) {
		string cmd1 = "plt.plot(";
		string s3 = arr_to_string_list(x, N1);

		string s4 = arr_to_string_list(z, N1);
		cmd1 += (s3 + "," + s4 + ", label='" + label2 + "')");
		PyRun_SimpleString(cmd1.c_str());
	}

	if (u != NULL) {
		string cmd2 = "plt.plot(";
		string s5 = arr_to_string_list(x, N1);

		string s6 = arr_to_string_list(u, N1);
		cmd2 += (s5 + "," + s6 + ", label='" + label3 + "')");
		PyRun_SimpleString(cmd2.c_str());
	}

	if (v != NULL) {
		string cmd3 = "plt.plot(";
		string s7 = arr_to_string_list(x, N1);

		string s8 = arr_to_string_list(v, N1);
		cmd3 += (s7 + "," + s8 + ", label='" + label4 + "')");
		PyRun_SimpleString(cmd3.c_str());
	}

	string cmd4 = "plt.title('";
	cmd4 += (title + "')");
	PyRun_SimpleString(cmd4.c_str());
	PyRun_SimpleString("plt.xlabel('time t')");
	PyRun_SimpleString("plt.legend()");


	PyRun_SimpleString("plt.show()");
}



void pythonInitial() {
	Py_Initialize(); /*初始化python解释器,告诉编译器要用的python编译器*/
	string path = ".";
	string chdir_cmd = string("sys.path.append(\"") + path + "\")";
	const char* cstr_cmd = chdir_cmd.c_str();
	PyRun_SimpleString("import sys");
	PyRun_SimpleString(cstr_cmd);
}

