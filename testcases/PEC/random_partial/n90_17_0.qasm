OPENQASM 2.0;
include "qelib1.inc";
qreg q[90];
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
cx q[19], q[35];
h q[43];
h q[88];
t q[30];
ccx q[3], q[16], q[36];
t q[42];
h q[71];
h q[67];
h q[79];
t q[39];
t q[73];
t q[59];
s q[58];
t q[81];
cx q[36], q[3];
ccx q[0], q[43], q[8];
h q[36];
s q[35];
ccx q[13], q[54], q[85];
ccx q[11], q[2], q[25];
cx q[8], q[67];
cx q[85], q[33];
cx q[89], q[44];
cx q[35], q[1];
t q[5];
ccx q[74], q[10], q[51];
h q[50];
s q[2];
h q[36];
h q[67];
t q[31];
h q[81];
ccx q[16], q[48], q[74];
cx q[34], q[49];
t q[84];
h q[46];
t q[58];
s q[16];
ccx q[82], q[35], q[4];
s q[42];
ccx q[77], q[45], q[58];
ccx q[58], q[7], q[79];
ccx q[51], q[85], q[66];
ccx q[36], q[31], q[48];
cx q[63], q[37];
cx q[44], q[25];
ccx q[29], q[45], q[14];
t q[63];
cx q[66], q[51];
s q[74];
cx q[71], q[67];
cx q[10], q[7];
cx q[24], q[7];
s q[70];
s q[9];
ccx q[49], q[3], q[16];
h q[69];
s q[41];
t q[2];
cx q[48], q[67];
cx q[88], q[7];
h q[52];
cx q[56], q[33];
ccx q[51], q[20], q[36];
ccx q[8], q[5], q[11];
h q[20];
t q[34];
t q[44];
t q[14];
h q[45];
s q[62];
s q[62];
s q[23];
cx q[23], q[39];
h q[23];
ccx q[30], q[24], q[46];
s q[85];
cx q[65], q[67];
t q[81];
s q[54];
cx q[26], q[28];
t q[56];
cx q[12], q[71];
s q[41];
s q[28];
h q[66];
s q[3];
cx q[84], q[6];
ccx q[40], q[34], q[65];
ccx q[39], q[47], q[68];
t q[46];
s q[56];
s q[31];
h q[82];
h q[72];
t q[60];
t q[64];
cx q[17], q[18];
ccx q[53], q[35], q[83];
s q[61];
h q[82];
ccx q[67], q[27], q[24];
t q[24];
s q[6];
t q[75];
h q[52];
h q[75];
cx q[0], q[88];
s q[78];
s q[77];
cx q[75], q[42];
t q[36];
t q[50];
s q[31];
h q[60];
h q[65];
h q[86];
t q[26];
t q[37];
cx q[22], q[26];
s q[50];
ccx q[62], q[26], q[2];
t q[21];
s q[44];
cx q[59], q[72];
t q[89];
cx q[71], q[64];
s q[24];
ccx q[45], q[54], q[80];
ccx q[0], q[54], q[82];
s q[89];
cx q[49], q[13];
s q[80];
t q[62];
cx q[27], q[56];
ccx q[7], q[60], q[83];
h q[70];
h q[70];
s q[12];
h q[18];
t q[88];
s q[69];
s q[68];
s q[47];
s q[76];
ccx q[38], q[68], q[73];
h q[70];
s q[26];
s q[69];
ccx q[36], q[35], q[26];
h q[67];
cx q[3], q[60];
h q[43];
t q[52];
cx q[45], q[80];
cx q[64], q[71];
h q[64];
ccx q[58], q[36], q[3];
t q[31];
t q[55];
h q[27];
h q[10];
s q[39];
s q[10];
cx q[73], q[18];
cx q[10], q[15];
t q[35];
cx q[20], q[46];
h q[5];
cx q[63], q[28];
t q[40];
h q[45];
s q[15];
cx q[43], q[12];
t q[3];
s q[69];
s q[62];
cx q[73], q[14];
cx q[83], q[35];
s q[33];
t q[20];
h q[81];
s q[42];
h q[45];
t q[36];
s q[80];
t q[20];
ccx q[38], q[53], q[61];
ccx q[74], q[61], q[88];
t q[73];
s q[3];
t q[69];
ccx q[65], q[73], q[57];
h q[43];
ccx q[20], q[38], q[70];
h q[33];
t q[46];
t q[28];
h q[51];
cx q[22], q[34];
s q[21];
ccx q[53], q[18], q[50];
s q[36];
cx q[60], q[2];
s q[11];
h q[55];
s q[63];
h q[72];
t q[34];
h q[84];
h q[83];
h q[73];
s q[38];
h q[73];
ccx q[34], q[50], q[25];
cx q[14], q[7];
h q[5];
t q[26];
cx q[82], q[16];
s q[58];
cx q[60], q[13];
t q[73];
s q[39];
cx q[11], q[61];
cx q[60], q[88];
h q[88];
h q[38];
cx q[33], q[68];
h q[25];
ccx q[30], q[3], q[89];
cx q[54], q[61];
s q[12];
cx q[34], q[0];
s q[87];
cx q[40], q[76];
cx q[88], q[73];
cx q[60], q[51];
ccx q[24], q[61], q[63];
h q[44];
s q[2];
t q[54];
s q[50];
t q[82];
ccx q[9], q[69], q[35];
cx q[66], q[13];
h q[61];
t q[14];
ccx q[31], q[38], q[28];
cx q[61], q[1];
t q[3];
t q[17];
ccx q[39], q[21], q[51];
cx q[82], q[7];
h q[13];
cx q[58], q[61];
cx q[29], q[49];
s q[65];
h q[81];
cx q[26], q[24];
cx q[45], q[85];
s q[33];
cx q[61], q[17];
ccx q[23], q[29], q[43];
ccx q[59], q[44], q[55];
ccx q[24], q[89], q[58];
t q[4];
t q[77];
cx q[50], q[70];
cx q[14], q[53];
h q[34];
x q[1];
cx q[2], q[3];
tdg q[4];
sdg q[4];
cx q[5], q[4];
sdg q[4];
z q[6];
z q[7];
cx q[9], q[8];
tdg q[8];
x q[8];
cx q[9], q[8];
x q[8];
x q[10];
cx q[10], q[11];
tdg q[11];
cx q[10], q[11];
x q[10];
cx q[13], q[14];
cx q[15], q[16];
tdg q[16];
sdg q[15];
cx q[15], q[16];
tdg q[18];
cx q[17], q[18];
sdg q[18];
tdg q[17];
cx q[20], q[19];
sdg q[19];
tdg q[19];
cx q[20], q[19];
tdg q[19];
sdg q[21];
cx q[21], q[22];
x q[23];
x q[24];
y q[26];
cx q[27], q[28];
cx q[30], q[29];
sdg q[29];
tdg q[33];
sdg q[33];
cx q[34], q[33];
tdg q[33];
cx q[35], q[36];
sdg q[38];
tdg q[37];
x q[39];
cx q[41], q[40];
sdg q[40];
tdg q[40];
cx q[41], q[40];
tdg q[40];
cx q[43], q[42];
tdg q[42];
ccx q[78], q[85], q[57];
cx q[49], q[87];
s q[89];
cx q[51], q[72];
h q[55];
s q[74];
ccx q[49], q[79], q[61];
cx q[84], q[67];
ccx q[64], q[62], q[59];
t q[62];
h q[82];
s q[84];
cx q[51], q[68];
h q[58];
ccx q[45], q[61], q[72];
cx q[60], q[48];
cx q[60], q[55];
s q[76];
h q[46];
t q[85];
h q[64];
ccx q[72], q[60], q[50];
ccx q[78], q[57], q[50];
ccx q[59], q[54], q[45];
t q[52];
h q[69];
t q[45];
t q[85];
cx q[87], q[48];
t q[47];
s q[60];
t q[53];
ccx q[62], q[81], q[78];
t q[53];
ccx q[79], q[69], q[45];
ccx q[77], q[85], q[84];
t q[48];
h q[83];
cx q[78], q[80];
cx q[48], q[54];
cx q[71], q[79];
ccx q[70], q[88], q[49];
t q[82];
t q[88];
ccx q[84], q[56], q[60];
