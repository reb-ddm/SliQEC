OPENQASM 2.0;
include "qelib1.inc";
qreg q[103];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
h q[30];
h q[31];
h q[32];
h q[33];
h q[34];
h q[35];
h q[36];
h q[37];
h q[38];
h q[39];
h q[40];
h q[41];
h q[42];
h q[43];
h q[44];
h q[45];
h q[46];
h q[47];
h q[48];
h q[49];
h q[50];
h q[51];
h q[52];
h q[53];
h q[54];
h q[55];
h q[56];
h q[57];
h q[58];
h q[59];
h q[60];
h q[61];
h q[62];
h q[63];
h q[64];
h q[65];
h q[66];
h q[67];
h q[68];
h q[69];
h q[70];
h q[71];
h q[72];
h q[73];
h q[74];
h q[75];
h q[76];
h q[77];
h q[78];
h q[79];
h q[80];
h q[81];
h q[82];
h q[83];
h q[84];
h q[85];
h q[86];
h q[87];
h q[88];
h q[89];
t q[20];
t q[22];
cx q[83], q[34];
s q[25];
h q[77];
ccx q[70], q[7], q[62];
s q[17];
t q[80];
h q[62];
t q[44];
t q[61];
ccx q[2], q[34], q[28];
s q[45];
cx q[27], q[26];
cx q[84], q[23];
ccx q[8], q[65], q[34];
cx q[0], q[45];
cx q[1], q[30];
cx q[5], q[11];
cx q[64], q[29];
ccx q[43], q[89], q[45];
cx q[65], q[7];
s q[81];
h q[67];
cx q[18], q[13];
cx q[48], q[30];
t q[0];
h q[33];
t q[14];
s q[21];
cx q[7], q[32];
s q[28];
s q[35];
h q[64];
ccx q[66], q[41], q[32];
cx q[12], q[46];
ccx q[37], q[88], q[65];
h q[51];
t q[49];
ccx q[65], q[61], q[33];
s q[53];
cx q[41], q[84];
h q[56];
cx q[2], q[85];
cx q[40], q[18];
ccx q[89], q[7], q[10];
ccx q[23], q[67], q[24];
h q[87];
ccx q[41], q[35], q[70];
t q[12];
s q[54];
ccx q[71], q[69], q[68];
s q[21];
s q[5];
cx q[84], q[14];
t q[79];
s q[32];
h q[48];
h q[80];
h q[15];
h q[48];
t q[65];
cx q[63], q[20];
cx q[65], q[66];
s q[44];
s q[36];
ccx q[81], q[9], q[64];
s q[26];
t q[89];
h q[30];
h q[82];
s q[63];
s q[21];
h q[55];
h q[51];
h q[84];
ccx q[74], q[36], q[55];
ccx q[34], q[22], q[74];
ccx q[58], q[3], q[84];
cx q[89], q[51];
s q[79];
h q[46];
t q[24];
cx q[57], q[38];
s q[20];
t q[54];
s q[38];
ccx q[0], q[66], q[63];
ccx q[13], q[33], q[54];
ccx q[15], q[74], q[5];
cx q[10], q[87];
s q[80];
t q[32];
t q[25];
s q[5];
cx q[44], q[15];
t q[8];
ccx q[68], q[5], q[15];
ccx q[46], q[64], q[13];
cx q[75], q[73];
ccx q[33], q[9], q[18];
cx q[33], q[66];
ccx q[68], q[76], q[57];
t q[77];
s q[56];
ccx q[30], q[11], q[12];
h q[55];
ccx q[1], q[5], q[86];
cx q[59], q[13];
h q[49];
s q[29];
ccx q[26], q[84], q[40];
ccx q[73], q[87], q[74];
h q[11];
ccx q[78], q[70], q[14];
h q[40];
s q[9];
ccx q[51], q[19], q[83];
cx q[47], q[31];
cx q[67], q[42];
cx q[36], q[57];
ccx q[36], q[87], q[25];
t q[45];
h q[42];
s q[59];
t q[28];
s q[47];
ccx q[17], q[55], q[53];
cx q[56], q[58];
t q[74];
h q[63];
ccx q[6], q[26], q[32];
ccx q[72], q[53], q[50];
ccx q[27], q[8], q[9];
cx q[3], q[59];
ccx q[75], q[10], q[2];
s q[24];
t q[35];
ccx q[22], q[34], q[16];
s q[14];
h q[9];
cx q[11], q[84];
cx q[55], q[60];
cx q[57], q[30];
cx q[33], q[66];
s q[75];
ccx q[33], q[0], q[45];
cx q[66], q[42];
ccx q[75], q[23], q[17];
ccx q[9], q[13], q[61];
t q[29];
cx q[63], q[38];
t q[31];
ccx q[50], q[34], q[60];
s q[48];
h q[33];
t q[59];
ccx q[56], q[69], q[70];
ccx q[40], q[5], q[59];
cx q[57], q[52];
cx q[40], q[34];
cx q[73], q[10];
s q[26];
cx q[13], q[16];
cx q[27], q[43];
h q[75];
s q[71];
t q[19];
cx q[55], q[56];
ccx q[66], q[36], q[40];
s q[84];
h q[54];
cx q[89], q[41];
s q[65];
t q[52];
s q[44];
ccx q[52], q[23], q[57];
cx q[44], q[86];
s q[71];
t q[62];
cx q[32], q[44];
cx q[39], q[18];
h q[79];
h q[30];
t q[67];
h q[4];
h q[1];
h q[5];
ccx q[84], q[76], q[83];
h q[12];
cx q[61], q[66];
t q[28];
ccx q[71], q[86], q[14];
h q[25];
ccx q[49], q[75], q[77];
cx q[82], q[27];
t q[21];
s q[71];
ccx q[38], q[61], q[79];
ccx q[1], q[8], q[3];
cx q[45], q[68];
h q[26];
ccx q[35], q[33], q[39];
s q[32];
t q[12];
h q[44];
h q[70];
h q[27];
h q[63];
s q[15];
cx q[27], q[36];
t q[25];
s q[87];
ccx q[83], q[45], q[62];
s q[59];
h q[31];
t q[28];
cx q[46], q[33];
h q[77];
cx q[81], q[15];
t q[31];
s q[6];
ccx q[15], q[47], q[79];
s q[22];
t q[25];
cx q[88], q[51];
s q[62];
h q[34];
t q[47];
ccx q[33], q[79], q[61];
t q[8];
cx q[35], q[50];
cx q[38], q[74];
s q[79];
h q[35];
h q[74];
ccx q[25], q[40], q[42];
h q[12];
ccx q[72], q[50], q[2];
h q[37];
h q[82];
t q[5];
ccx q[14], q[67], q[89];
t q[77];
ccx q[50], q[8], q[47];
cx q[4], q[73];
s q[78];
t q[24];
t q[45];
ccx q[86], q[7], q[19];
ccx q[68], q[4], q[81];
ccx q[50], q[24], q[51];
s q[48];
t q[89];
h q[62];
ccx q[40], q[68], q[70];
h q[82];
t q[20];
s q[17];
cx q[58], q[28];
s q[23];
t q[64];
s q[51];
t q[67];
ccx q[16], q[49], q[35];
s q[75];
h q[71];
ccx q[24], q[25], q[50];
ccx q[70], q[41], q[40];
ccx q[39], q[49], q[22];
z q[0];
tdg q[1];
cx q[1], q[2];
tdg q[4];
sdg q[5];
cx q[4], q[5];
tdg q[4];
sdg q[6];
cx q[9], q[8];
sdg q[8];
z q[10];
x q[11];
x q[12];
cx q[13], q[14];
tdg q[14];
cx q[13], q[14];
tdg q[13];
sdg q[13];
sdg q[15];
sdg q[17];
cx q[18], q[17];
tdg q[17];
sdg q[19];
tdg q[21];
tdg q[23];
cx q[23], q[24];
tdg q[26];
cx q[25], q[26];
tdg q[25];
z q[27];
sdg q[29];
cx q[28], q[29];
sdg q[29];
tdg q[28];
cx q[30], q[31];
tdg q[31];
sdg q[30];
cx q[30], q[31];
cx q[32], q[33];
cx q[35], q[34];
tdg q[34];
cx q[36], q[37];
tdg q[36];
sdg q[39];
cx q[41], q[42];
cx q[43], q[44];
sdg q[44];
sdg q[43];
h q[49];
cx q[55], q[63];
ccx q[65], q[79], q[49];
ccx q[70], q[45], q[61];
ccx q[87], q[60], q[59];
ccx q[55], q[69], q[65];
ccx q[74], q[85], q[82];
s q[84];
t q[66];
h q[59];
t q[66];
h q[51];
s q[82];
ccx q[47], q[66], q[83];
ccx q[46], q[71], q[56];
ccx q[88], q[59], q[49];
t q[70];
t q[87];
cx q[56], q[57];
h q[57];
cx q[89], q[69];
s q[79];
t q[85];
t q[65];
s q[56];
cx q[76], q[54];
cx q[87], q[74];
s q[89];
t q[76];
t q[54];
ccx q[82], q[66], q[73];
cx q[81], q[59];
ccx q[52], q[45], q[50];
ccx q[45], q[63], q[70];
ccx q[71], q[80], q[64];
s q[78];
cx q[74], q[62];
s q[72];
s q[70];
t q[62];
cx q[60], q[89];
ccx q[70], q[45], q[52];
t q[78];
t q[89];
h q[45];
cx q[90], q[88];
cx q[91], q[13];
cx q[92], q[22];
cx q[93], q[31];
cx q[94], q[10];
cx q[95], q[74];
cx q[96], q[19];
cx q[97], q[26];
cx q[98], q[28];
cx q[99], q[39];
cx q[100], q[18];
cx q[101], q[7];
cx q[102], q[29];