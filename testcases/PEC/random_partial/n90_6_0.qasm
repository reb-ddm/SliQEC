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
t q[6];
h q[37];
h q[62];
s q[74];
cx q[65], q[74];
cx q[11], q[0];
h q[42];
t q[29];
cx q[62], q[74];
h q[83];
s q[59];
cx q[24], q[65];
ccx q[33], q[40], q[83];
t q[34];
cx q[73], q[87];
cx q[10], q[65];
t q[82];
s q[2];
cx q[86], q[76];
s q[74];
ccx q[44], q[22], q[7];
ccx q[49], q[61], q[75];
cx q[29], q[23];
h q[72];
cx q[30], q[77];
t q[13];
ccx q[6], q[13], q[24];
t q[69];
h q[70];
cx q[33], q[5];
t q[83];
t q[57];
cx q[5], q[89];
ccx q[24], q[49], q[82];
ccx q[5], q[82], q[75];
ccx q[32], q[35], q[69];
t q[50];
h q[60];
h q[52];
h q[13];
t q[32];
cx q[81], q[33];
cx q[54], q[47];
ccx q[64], q[38], q[5];
s q[19];
cx q[67], q[56];
ccx q[83], q[27], q[74];
t q[16];
s q[48];
h q[34];
h q[58];
ccx q[66], q[50], q[77];
t q[54];
t q[61];
h q[5];
h q[32];
t q[58];
cx q[23], q[70];
ccx q[33], q[48], q[87];
h q[55];
h q[21];
s q[72];
s q[5];
cx q[10], q[34];
h q[7];
h q[36];
t q[44];
h q[83];
s q[80];
s q[20];
s q[79];
ccx q[6], q[81], q[37];
ccx q[74], q[10], q[33];
s q[53];
h q[74];
t q[14];
h q[45];
h q[28];
ccx q[17], q[74], q[82];
t q[23];
ccx q[66], q[35], q[55];
cx q[71], q[57];
cx q[31], q[89];
s q[78];
h q[61];
s q[71];
h q[56];
t q[53];
t q[72];
s q[70];
t q[46];
t q[2];
cx q[13], q[4];
h q[70];
t q[23];
h q[63];
ccx q[73], q[79], q[0];
cx q[45], q[31];
s q[16];
h q[83];
ccx q[83], q[16], q[23];
s q[37];
h q[54];
cx q[51], q[1];
t q[66];
ccx q[77], q[37], q[1];
t q[17];
ccx q[5], q[42], q[76];
t q[73];
t q[35];
t q[1];
t q[50];
ccx q[59], q[19], q[10];
ccx q[20], q[29], q[34];
t q[85];
s q[60];
t q[10];
cx q[52], q[53];
cx q[71], q[65];
s q[68];
ccx q[69], q[81], q[34];
h q[69];
t q[33];
t q[39];
s q[29];
h q[55];
ccx q[70], q[79], q[59];
cx q[57], q[34];
t q[64];
h q[56];
h q[86];
ccx q[34], q[5], q[80];
s q[11];
s q[12];
cx q[14], q[35];
t q[51];
t q[86];
t q[10];
cx q[43], q[48];
cx q[59], q[78];
s q[26];
s q[35];
cx q[71], q[81];
ccx q[43], q[51], q[36];
cx q[68], q[34];
t q[85];
ccx q[7], q[24], q[15];
h q[29];
t q[33];
t q[61];
ccx q[62], q[52], q[76];
s q[57];
t q[35];
cx q[39], q[51];
ccx q[85], q[13], q[59];
ccx q[9], q[55], q[39];
cx q[51], q[66];
cx q[26], q[45];
cx q[4], q[16];
t q[4];
ccx q[73], q[4], q[39];
cx q[44], q[32];
cx q[34], q[72];
ccx q[79], q[29], q[89];
cx q[25], q[7];
h q[31];
ccx q[17], q[24], q[43];
s q[86];
h q[51];
h q[11];
ccx q[30], q[62], q[21];
ccx q[55], q[31], q[24];
s q[61];
ccx q[10], q[66], q[40];
ccx q[65], q[0], q[2];
t q[77];
cx q[52], q[51];
s q[42];
ccx q[83], q[33], q[4];
s q[47];
ccx q[83], q[3], q[30];
h q[7];
s q[14];
s q[6];
t q[36];
cx q[58], q[87];
s q[42];
ccx q[73], q[47], q[20];
t q[89];
h q[11];
s q[73];
ccx q[19], q[16], q[58];
cx q[57], q[15];
cx q[2], q[17];
t q[60];
s q[14];
t q[36];
ccx q[4], q[87], q[69];
t q[65];
t q[41];
cx q[58], q[39];
t q[75];
ccx q[29], q[34], q[69];
s q[16];
s q[79];
h q[18];
s q[67];
h q[44];
h q[7];
ccx q[27], q[67], q[21];
s q[23];
h q[10];
ccx q[21], q[3], q[14];
t q[45];
h q[81];
cx q[47], q[22];
ccx q[83], q[69], q[58];
h q[54];
ccx q[24], q[41], q[47];
h q[38];
t q[5];
ccx q[69], q[35], q[1];
cx q[2], q[14];
ccx q[4], q[84], q[55];
t q[56];
s q[83];
ccx q[34], q[37], q[86];
t q[52];
h q[35];
s q[55];
cx q[8], q[21];
ccx q[35], q[56], q[5];
ccx q[41], q[36], q[34];
h q[62];
ccx q[62], q[37], q[21];
t q[1];
s q[58];
t q[22];
ccx q[59], q[17], q[80];
t q[64];
cx q[15], q[2];
cx q[81], q[7];
s q[60];
cx q[33], q[76];
t q[75];
s q[75];
ccx q[35], q[13], q[17];
t q[34];
s q[14];
s q[60];
h q[45];
t q[20];
ccx q[49], q[22], q[66];
ccx q[27], q[32], q[55];
t q[66];
s q[62];
t q[1];
s q[22];
cx q[4], q[20];
ccx q[36], q[54], q[23];
ccx q[12], q[17], q[39];
ccx q[45], q[21], q[9];
cx q[46], q[75];
cx q[72], q[68];
ccx q[47], q[13], q[29];
cx q[54], q[81];
cx q[5], q[28];
h q[25];
h q[83];
h q[46];
cx q[1], q[2];
tdg q[1];
tdg q[4];
tdg q[3];
cx q[3], q[4];
tdg q[3];
z q[5];
y q[6];
sdg q[7];
cx q[7], q[8];
cx q[10], q[9];
sdg q[9];
cx q[10], q[9];
sdg q[9];
x q[11];
x q[12];
sdg q[13];
cx q[15], q[16];
tdg q[17];
tdg q[19];
tdg q[22];
tdg q[21];
cx q[21], q[22];
z q[23];
sdg q[24];
cx q[26], q[27];
tdg q[26];
tdg q[29];
cx q[28], q[29];
tdg q[28];
sdg q[30];
cx q[31], q[30];
tdg q[30];
sdg q[30];
z q[32];
cx q[35], q[34];
sdg q[34];
sdg q[35];
cx q[35], q[34];
sdg q[34];
cx q[36], q[37];
tdg q[36];
y q[38];
tdg q[39];
sdg q[39];
z q[41];
tdg q[43];
sdg q[42];
cx q[43], q[42];
tdg q[42];
h q[79];
cx q[48], q[78];
h q[70];
t q[72];
t q[68];
h q[87];
t q[54];
t q[80];
t q[58];
t q[59];
t q[81];
t q[48];
cx q[74], q[64];
h q[48];
h q[56];
cx q[85], q[68];
t q[58];
h q[78];
h q[80];
ccx q[45], q[74], q[66];
t q[75];
ccx q[51], q[89], q[63];
cx q[49], q[57];
cx q[55], q[70];
ccx q[61], q[65], q[63];
h q[52];
t q[74];
s q[69];
t q[83];
ccx q[74], q[55], q[73];
cx q[51], q[80];
cx q[49], q[45];
s q[84];
t q[53];
ccx q[66], q[60], q[73];
h q[85];
h q[50];
h q[49];
h q[88];
s q[59];
ccx q[77], q[56], q[66];
s q[48];
t q[86];
s q[68];
cx q[51], q[59];
