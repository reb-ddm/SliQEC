OPENQASM 2.0;
include "qelib1.inc";
qreg q[89];
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
t q[22];
h q[36];
cx q[38], q[39];
ccx q[22], q[65], q[18];
h q[55];
t q[16];
s q[56];
cx q[34], q[0];
ccx q[7], q[33], q[59];
ccx q[59], q[68], q[39];
ccx q[15], q[48], q[46];
t q[59];
ccx q[20], q[65], q[48];
s q[74];
ccx q[3], q[10], q[77];
cx q[20], q[2];
h q[7];
ccx q[38], q[17], q[41];
ccx q[58], q[56], q[2];
h q[59];
h q[45];
h q[71];
cx q[3], q[5];
cx q[3], q[28];
cx q[59], q[78];
cx q[46], q[38];
ccx q[50], q[52], q[20];
s q[14];
ccx q[24], q[21], q[60];
t q[35];
cx q[25], q[56];
cx q[23], q[42];
s q[79];
s q[64];
cx q[77], q[62];
s q[76];
cx q[66], q[73];
cx q[32], q[42];
cx q[71], q[28];
s q[34];
t q[41];
s q[15];
cx q[34], q[20];
cx q[6], q[50];
cx q[51], q[28];
cx q[35], q[52];
t q[79];
h q[27];
t q[6];
ccx q[71], q[67], q[12];
h q[50];
cx q[9], q[15];
ccx q[70], q[51], q[57];
s q[32];
h q[43];
ccx q[49], q[6], q[13];
s q[75];
cx q[47], q[43];
s q[55];
cx q[25], q[18];
s q[18];
s q[54];
ccx q[11], q[54], q[19];
h q[42];
t q[48];
h q[74];
cx q[70], q[37];
cx q[46], q[50];
cx q[3], q[57];
t q[78];
ccx q[77], q[32], q[38];
h q[63];
s q[67];
cx q[21], q[42];
s q[24];
s q[27];
t q[57];
t q[20];
h q[13];
h q[8];
cx q[31], q[75];
t q[57];
cx q[62], q[33];
cx q[43], q[17];
cx q[66], q[57];
cx q[62], q[61];
t q[73];
ccx q[23], q[57], q[77];
s q[21];
h q[33];
cx q[50], q[22];
s q[18];
ccx q[61], q[34], q[2];
ccx q[10], q[19], q[73];
h q[45];
t q[77];
ccx q[35], q[27], q[55];
s q[49];
t q[37];
t q[74];
s q[13];
ccx q[47], q[16], q[42];
h q[0];
s q[34];
ccx q[20], q[78], q[57];
t q[76];
ccx q[16], q[56], q[46];
h q[47];
t q[47];
cx q[78], q[63];
t q[54];
t q[12];
t q[54];
h q[0];
t q[73];
ccx q[11], q[23], q[73];
s q[6];
s q[4];
h q[38];
t q[77];
cx q[77], q[34];
ccx q[30], q[3], q[60];
cx q[54], q[38];
t q[57];
ccx q[50], q[36], q[34];
h q[63];
ccx q[36], q[2], q[49];
ccx q[64], q[6], q[34];
s q[47];
ccx q[6], q[29], q[54];
ccx q[42], q[31], q[56];
cx q[36], q[69];
h q[56];
ccx q[42], q[53], q[7];
t q[78];
s q[69];
s q[60];
ccx q[74], q[33], q[30];
s q[22];
ccx q[32], q[79], q[35];
ccx q[54], q[35], q[51];
ccx q[3], q[66], q[43];
t q[10];
cx q[56], q[73];
h q[7];
t q[15];
cx q[55], q[26];
cx q[73], q[31];
cx q[10], q[14];
ccx q[48], q[50], q[9];
cx q[74], q[24];
ccx q[27], q[35], q[3];
ccx q[38], q[22], q[76];
t q[22];
s q[50];
t q[64];
s q[14];
cx q[1], q[41];
ccx q[14], q[27], q[20];
t q[10];
ccx q[56], q[0], q[68];
cx q[48], q[56];
cx q[4], q[12];
t q[20];
s q[27];
ccx q[41], q[63], q[17];
cx q[25], q[68];
t q[42];
cx q[66], q[63];
cx q[13], q[78];
cx q[37], q[17];
h q[51];
h q[50];
s q[78];
cx q[6], q[21];
s q[32];
cx q[55], q[4];
s q[21];
t q[18];
s q[2];
cx q[51], q[20];
ccx q[5], q[73], q[59];
h q[14];
ccx q[77], q[3], q[23];
ccx q[0], q[24], q[12];
s q[77];
s q[46];
cx q[39], q[61];
s q[8];
ccx q[27], q[53], q[33];
ccx q[0], q[13], q[26];
s q[44];
s q[2];
s q[77];
s q[45];
s q[70];
ccx q[43], q[26], q[77];
t q[63];
cx q[1], q[77];
t q[68];
ccx q[21], q[18], q[59];
t q[49];
ccx q[23], q[14], q[76];
ccx q[78], q[32], q[70];
t q[54];
t q[11];
t q[57];
ccx q[8], q[31], q[38];
ccx q[8], q[18], q[7];
cx q[53], q[22];
h q[63];
ccx q[44], q[39], q[24];
s q[28];
cx q[76], q[17];
ccx q[1], q[78], q[26];
s q[64];
cx q[42], q[34];
ccx q[62], q[25], q[68];
t q[43];
cx q[53], q[51];
ccx q[2], q[24], q[46];
ccx q[41], q[58], q[47];
ccx q[27], q[6], q[25];
s q[57];
t q[9];
ccx q[17], q[71], q[29];
t q[7];
h q[42];
ccx q[33], q[41], q[3];
cx q[72], q[6];
s q[51];
s q[54];
t q[12];
s q[41];
t q[65];
h q[77];
s q[9];
t q[13];
ccx q[1], q[68], q[3];
cx q[0], q[26];
sdg q[0];
tdg q[1];
tdg q[0];
cx q[0], q[1];
cx q[2], q[3];
tdg q[3];
cx q[2], q[3];
sdg q[2];
cx q[4], q[5];
tdg q[4];
sdg q[6];
cx q[6], q[7];
y q[8];
cx q[9], q[10];
tdg q[10];
cx q[9], q[10];
sdg q[9];
cx q[11], q[12];
tdg q[12];
tdg q[11];
sdg q[11];
cx q[11], q[12];
cx q[13], q[14];
tdg q[13];
sdg q[14];
cx q[13], q[14];
sdg q[13];
sdg q[15];
cx q[16], q[15];
sdg q[15];
tdg q[19];
tdg q[18];
cx q[19], q[18];
sdg q[18];
z q[20];
cx q[21], q[22];
tdg q[21];
cx q[24], q[23];
tdg q[24];
tdg q[23];
cx q[24], q[23];
tdg q[23];
tdg q[26];
cx q[25], q[26];
sdg q[26];
tdg q[25];
tdg q[27];
tdg q[30];
cx q[29], q[30];
tdg q[29];
x q[31];
cx q[33], q[32];
sdg q[33];
tdg q[32];
cx q[35], q[34];
tdg q[34];
sdg q[36];
tdg q[36];
cx q[37], q[36];
tdg q[36];
z q[38];
y q[39];
cx q[51], q[73];
ccx q[40], q[52], q[65];
ccx q[76], q[70], q[48];
t q[57];
ccx q[52], q[55], q[72];
s q[56];
s q[45];
s q[72];
s q[42];
cx q[40], q[64];
cx q[78], q[77];
h q[47];
ccx q[48], q[62], q[43];
ccx q[61], q[69], q[78];
h q[62];
cx q[78], q[67];
t q[51];
s q[58];
ccx q[40], q[56], q[58];
cx q[72], q[71];
s q[40];
cx q[41], q[63];
s q[48];
ccx q[56], q[52], q[48];
ccx q[60], q[64], q[57];
s q[45];
t q[58];
s q[58];
h q[47];
ccx q[56], q[65], q[78];
h q[52];
cx q[78], q[66];
t q[77];
cx q[75], q[67];
cx q[54], q[71];
s q[58];
ccx q[77], q[44], q[63];
t q[56];
t q[62];
h q[53];
cx q[80], q[64];
cx q[81], q[27];
cx q[82], q[13];
cx q[83], q[77];
cx q[84], q[8];
cx q[85], q[41];
cx q[86], q[10];
cx q[87], q[44];
cx q[88], q[5];
