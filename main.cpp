#include "EMG.h"
extern "C" { FILE __iob_func[3] = { *stdin,*stdout,*stderr }; }
int main() {
	
	mat output,result_map;
	output.load("output.csv", csv_ascii);
	//output.print();

	result_map = EMG(output);
	result_map.save("result_map.csv", csv_ascii);


	return 0;
}