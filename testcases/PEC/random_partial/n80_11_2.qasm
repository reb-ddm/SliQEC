OPENQASM 2.0;
include "qelib1.inc";
qreg q[92];
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
s q[34];
s q[21];
ccx q[76], q[9], q[19];
t q[2];
cx q[6], q[33];
cx q[15], q[47];
s q[39];
cx q[39], q[32];
t q[39];
cx q[73], q[12];
h q[72];
cx q[21], q[75];
t q[61];
s q[55];
s q[47];
t q[62];
s q[64];
h q[47];
t q[17];
cx q[18], q[48];
cx q[32], q[29];
cx q[7], q[78];
cx q[4], q[14];
cx q[16], q[8];
t q[3];
t q[56];
s q[70];
s q[57];
h q[30];
cx q[30], q[26];
ccx q[35], q[48], q[65];
s q[52];
t q[39];
h q[4];
s q[34];
ccx q[38], q[6], q[7];
t q[15];
s q[25];
s q[21];
s q[57];
cx q[31], q[44];
cx q[50], q[38];
h q[65];
t q[22];
cx q[57], q[26];
h q[76];
ccx q[27], q[37], q[64];
t q[48];
h q[18];
cx q[53], q[33];
ccx q[28], q[21], q[5];
h q[63];
cx q[69], q[3];
h q[5];
ccx q[7], q[30], q[38];
s q[6];
ccx q[68], q[38], q[31];
s q[23];
ccx q[34], q[37], q[43];
t q[65];
ccx q[32], q[1], q[12];
s q[70];
h q[33];
s q[9];
cx q[21], q[73];
ccx q[56], q[68], q[37];
cx q[4], q[71];
h q[61];
h q[15];
h q[24];
cx q[45], q[38];
t q[15];
s q[11];
t q[33];
ccx q[20], q[41], q[60];
h q[43];
cx q[68], q[3];
h q[6];
cx q[62], q[21];
s q[11];
h q[4];
t q[37];
t q[5];
cx q[76], q[50];
ccx q[26], q[16], q[57];
t q[53];
t q[18];
cx q[35], q[59];
cx q[7], q[51];
cx q[55], q[32];
ccx q[27], q[45], q[67];
s q[72];
cx q[48], q[14];
s q[56];
cx q[24], q[66];
s q[67];
t q[40];
s q[37];
h q[75];
ccx q[77], q[58], q[68];
ccx q[0], q[77], q[73];
t q[35];
h q[3];
h q[5];
cx q[30], q[61];
cx q[76], q[77];
ccx q[45], q[54], q[13];
ccx q[54], q[27], q[49];
t q[6];
cx q[1], q[61];
s q[1];
s q[27];
cx q[54], q[6];
h q[9];
t q[54];
cx q[50], q[67];
ccx q[41], q[51], q[24];
ccx q[57], q[50], q[77];
cx q[22], q[71];
t q[32];
s q[30];
s q[55];
cx q[7], q[27];
t q[7];
s q[38];
ccx q[33], q[77], q[17];
s q[61];
s q[15];
t q[8];
t q[62];
t q[59];
t q[79];
ccx q[35], q[79], q[2];
ccx q[38], q[55], q[16];
t q[56];
ccx q[57], q[68], q[12];
h q[35];
h q[35];
cx q[28], q[10];
s q[54];
t q[63];
ccx q[28], q[49], q[4];
t q[40];
s q[61];
cx q[63], q[11];
h q[75];
ccx q[19], q[44], q[36];
ccx q[44], q[20], q[40];
s q[57];
t q[5];
ccx q[15], q[42], q[37];
s q[78];
s q[32];
h q[76];
h q[67];
t q[73];
t q[54];
h q[29];
s q[8];
ccx q[48], q[23], q[3];
h q[68];
h q[29];
s q[40];
t q[42];
cx q[66], q[18];
cx q[35], q[74];
cx q[16], q[4];
s q[65];
s q[56];
s q[51];
h q[74];
s q[8];
ccx q[18], q[0], q[16];
t q[3];
s q[26];
t q[9];
h q[15];
s q[26];
s q[61];
t q[23];
t q[40];
h q[37];
cx q[14], q[10];
t q[75];
cx q[46], q[79];
cx q[40], q[1];
ccx q[54], q[14], q[58];
h q[71];
t q[72];
s q[51];
cx q[6], q[22];
h q[25];
h q[56];
ccx q[24], q[45], q[16];
t q[41];
s q[75];
s q[62];
h q[39];
ccx q[66], q[41], q[26];
s q[47];
cx q[71], q[24];
s q[59];
t q[69];
h q[12];
t q[20];
h q[1];
t q[33];
t q[30];
t q[65];
cx q[65], q[42];
t q[70];
h q[12];
cx q[22], q[1];
t q[42];
ccx q[8], q[76], q[58];
t q[46];
cx q[61], q[38];
h q[14];
h q[51];
s q[34];
s q[76];
ccx q[12], q[41], q[62];
cx q[51], q[71];
ccx q[4], q[45], q[37];
ccx q[38], q[69], q[72];
t q[10];
ccx q[69], q[15], q[12];
t q[62];
cx q[47], q[62];
t q[29];
h q[77];
h q[2];
t q[67];
t q[28];
t q[2];
s q[1];
t q[31];
t q[68];
h q[52];
h q[61];
tdg q[2];
cx q[1], q[2];
sdg q[1];
tdg q[3];
cx q[3], q[4];
tdg q[4];
cx q[3], q[4];
sdg q[3];
tdg q[6];
tdg q[7];
cx q[7], q[6];
tdg q[6];
x q[8];
sdg q[10];
cx q[10], q[11];
tdg q[12];
sdg q[15];
cx q[14], q[15];
tdg q[14];
tdg q[16];
tdg q[18];
x q[20];
tdg q[23];
sdg q[22];
tdg q[22];
cx q[22], q[23];
cx q[26], q[25];
tdg q[25];
sdg q[26];
cx q[26], q[25];
tdg q[25];
cx q[28], q[27];
tdg q[27];
z q[29];
tdg q[31];
cx q[30], q[31];
tdg q[30];
sdg q[32];
cx q[32], q[33];
tdg q[33];
tdg q[32];
cx q[32], q[33];
y q[34];
sdg q[36];
cx q[36], q[35];
tdg q[35];
cx q[38], q[37];
sdg q[37];
cx q[43], q[73];
cx q[70], q[64];
s q[60];
cx q[53], q[76];
h q[68];
cx q[61], q[73];
s q[45];
cx q[54], q[55];
t q[44];
h q[66];
s q[48];
cx q[69], q[47];
cx q[70], q[46];
h q[63];
s q[64];
cx q[59], q[78];
ccx q[51], q[57], q[55];
h q[53];
cx q[41], q[61];
ccx q[63], q[53], q[71];
h q[49];
t q[43];
t q[58];
ccx q[47], q[72], q[48];
ccx q[54], q[44], q[51];
cx q[41], q[52];
t q[45];
t q[43];
h q[66];
t q[69];
t q[65];
cx q[46], q[79];
cx q[44], q[66];
s q[68];
h q[69];
h q[51];
cx q[57], q[49];
cx q[78], q[69];
h q[58];
s q[45];
cx q[80], q[32];
cx q[81], q[55];
cx q[82], q[78];
cx q[83], q[35];
cx q[84], q[76];
cx q[85], q[26];
cx q[86], q[72];
cx q[87], q[28];
cx q[88], q[24];
cx q[89], q[14];
cx q[90], q[74];
cx q[91], q[50];
