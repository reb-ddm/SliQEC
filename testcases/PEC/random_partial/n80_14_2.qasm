OPENQASM 2.0;
include "qelib1.inc";
qreg q[88];
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
cx q[45], q[56];
cx q[20], q[15];
h q[20];
s q[73];
s q[45];
t q[28];
s q[9];
cx q[33], q[71];
s q[70];
cx q[41], q[77];
s q[60];
h q[20];
t q[37];
cx q[76], q[74];
s q[62];
s q[53];
cx q[75], q[5];
h q[59];
t q[57];
h q[21];
s q[55];
t q[21];
h q[62];
cx q[21], q[10];
ccx q[54], q[49], q[65];
cx q[25], q[11];
t q[36];
ccx q[31], q[42], q[43];
ccx q[15], q[63], q[64];
ccx q[77], q[72], q[37];
h q[8];
cx q[1], q[63];
ccx q[27], q[0], q[40];
ccx q[73], q[42], q[45];
s q[78];
s q[79];
ccx q[59], q[5], q[7];
h q[66];
cx q[28], q[19];
h q[16];
ccx q[24], q[21], q[28];
s q[40];
s q[66];
ccx q[10], q[7], q[29];
s q[54];
t q[4];
t q[61];
h q[66];
t q[12];
h q[14];
ccx q[25], q[45], q[22];
cx q[29], q[28];
ccx q[52], q[25], q[2];
h q[70];
cx q[46], q[60];
t q[52];
cx q[47], q[27];
ccx q[32], q[78], q[51];
cx q[58], q[39];
h q[52];
h q[6];
ccx q[32], q[70], q[37];
h q[34];
s q[31];
t q[49];
h q[60];
t q[72];
cx q[52], q[79];
cx q[56], q[9];
s q[51];
h q[26];
cx q[22], q[34];
cx q[15], q[6];
s q[9];
t q[49];
t q[47];
ccx q[69], q[24], q[66];
s q[19];
t q[22];
ccx q[9], q[50], q[31];
ccx q[44], q[55], q[46];
t q[36];
cx q[13], q[74];
t q[66];
ccx q[53], q[10], q[22];
cx q[16], q[56];
ccx q[58], q[9], q[6];
h q[35];
h q[26];
h q[41];
ccx q[61], q[52], q[74];
s q[22];
s q[34];
ccx q[59], q[51], q[7];
ccx q[14], q[48], q[67];
ccx q[15], q[14], q[21];
h q[54];
s q[47];
h q[71];
cx q[76], q[65];
cx q[7], q[31];
h q[39];
ccx q[0], q[24], q[76];
ccx q[16], q[33], q[66];
t q[62];
s q[15];
t q[16];
cx q[12], q[56];
s q[27];
s q[45];
t q[10];
s q[8];
cx q[34], q[22];
cx q[52], q[48];
ccx q[39], q[5], q[20];
ccx q[34], q[39], q[70];
ccx q[9], q[61], q[33];
h q[20];
h q[76];
cx q[46], q[54];
t q[15];
t q[0];
h q[53];
ccx q[6], q[21], q[15];
t q[4];
cx q[24], q[11];
cx q[63], q[14];
t q[5];
ccx q[37], q[10], q[1];
cx q[61], q[20];
s q[75];
t q[20];
cx q[41], q[7];
s q[46];
t q[37];
cx q[63], q[54];
h q[77];
cx q[68], q[41];
t q[74];
h q[9];
ccx q[56], q[16], q[25];
h q[72];
t q[17];
s q[69];
h q[28];
cx q[54], q[28];
ccx q[26], q[78], q[72];
t q[18];
ccx q[31], q[2], q[72];
s q[68];
t q[46];
t q[64];
cx q[25], q[20];
cx q[21], q[29];
t q[57];
ccx q[9], q[71], q[42];
ccx q[64], q[8], q[73];
cx q[61], q[34];
s q[68];
h q[52];
h q[54];
ccx q[36], q[27], q[51];
s q[7];
t q[21];
s q[15];
t q[73];
t q[3];
ccx q[2], q[9], q[6];
cx q[30], q[68];
s q[12];
ccx q[40], q[23], q[9];
s q[64];
cx q[76], q[18];
cx q[31], q[29];
t q[68];
t q[64];
ccx q[57], q[31], q[16];
cx q[21], q[63];
s q[13];
s q[70];
ccx q[64], q[58], q[52];
h q[17];
ccx q[59], q[47], q[77];
s q[77];
t q[67];
h q[37];
cx q[49], q[45];
cx q[7], q[51];
h q[27];
t q[18];
ccx q[54], q[32], q[43];
ccx q[37], q[21], q[6];
ccx q[63], q[43], q[13];
h q[1];
cx q[44], q[63];
t q[61];
ccx q[40], q[77], q[18];
ccx q[72], q[75], q[67];
h q[23];
ccx q[20], q[19], q[66];
s q[36];
cx q[67], q[25];
h q[54];
t q[58];
t q[29];
s q[37];
h q[11];
t q[78];
t q[41];
h q[3];
s q[29];
t q[78];
h q[69];
t q[41];
s q[73];
ccx q[67], q[20], q[15];
ccx q[36], q[35], q[21];
s q[58];
cx q[0], q[10];
cx q[57], q[74];
s q[6];
t q[43];
s q[3];
cx q[55], q[26];
ccx q[24], q[4], q[28];
s q[71];
h q[61];
cx q[64], q[1];
ccx q[9], q[62], q[55];
h q[77];
t q[31];
h q[44];
cx q[13], q[41];
ccx q[14], q[35], q[18];
h q[69];
h q[42];
cx q[77], q[24];
h q[5];
h q[30];
s q[61];
cx q[1], q[0];
tdg q[0];
x q[2];
z q[3];
sdg q[4];
sdg q[7];
tdg q[6];
z q[8];
cx q[9], q[10];
tdg q[10];
cx q[9], q[10];
sdg q[10];
tdg q[9];
cx q[11], q[12];
tdg q[12];
tdg q[11];
cx q[11], q[12];
sdg q[11];
sdg q[13];
cx q[14], q[13];
tdg q[14];
tdg q[13];
sdg q[16];
cx q[15], q[16];
tdg q[16];
sdg q[15];
cx q[15], q[16];
cx q[19], q[18];
tdg q[19];
sdg q[18];
z q[20];
cx q[21], q[22];
sdg q[21];
cx q[23], q[24];
tdg q[24];
tdg q[23];
cx q[23], q[24];
tdg q[23];
sdg q[26];
cx q[25], q[26];
sdg q[26];
tdg q[25];
cx q[28], q[27];
sdg q[27];
x q[29];
cx q[29], q[30];
sdg q[30];
cx q[29], q[30];
x q[29];
sdg q[32];
sdg q[31];
cx q[31], q[32];
y q[33];
cx q[35], q[34];
tdg q[34];
cx q[35], q[34];
tdg q[34];
y q[36];
z q[37];
y q[38];
ccx q[78], q[77], q[44];
cx q[77], q[58];
s q[50];
ccx q[57], q[69], q[40];
s q[46];
t q[46];
ccx q[57], q[79], q[54];
cx q[40], q[54];
ccx q[65], q[68], q[42];
t q[76];
cx q[64], q[61];
s q[45];
s q[53];
h q[77];
ccx q[70], q[74], q[52];
h q[44];
t q[52];
ccx q[43], q[65], q[73];
ccx q[53], q[71], q[79];
t q[56];
h q[53];
h q[43];
t q[60];
s q[51];
s q[40];
t q[61];
t q[47];
t q[44];
cx q[42], q[54];
s q[40];
h q[70];
ccx q[45], q[43], q[42];
ccx q[45], q[64], q[53];
h q[64];
cx q[50], q[66];
ccx q[55], q[74], q[76];
h q[47];
t q[79];
h q[47];
s q[46];
cx q[80], q[19];
cx q[81], q[55];
cx q[82], q[26];
cx q[83], q[7];
cx q[84], q[18];
cx q[85], q[34];
cx q[86], q[67];
cx q[87], q[64];