OPENQASM 2.0;
include "qelib1.inc";
qreg q[105];
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
cx q[20], q[28];
cx q[51], q[58];
ccx q[15], q[35], q[88];
t q[80];
cx q[31], q[65];
s q[52];
cx q[38], q[66];
ccx q[55], q[5], q[39];
cx q[57], q[21];
t q[57];
t q[7];
s q[58];
cx q[49], q[82];
s q[82];
s q[72];
ccx q[36], q[57], q[65];
h q[55];
s q[57];
cx q[55], q[0];
ccx q[30], q[2], q[23];
s q[75];
cx q[1], q[16];
t q[25];
cx q[77], q[16];
ccx q[19], q[3], q[54];
ccx q[62], q[71], q[14];
cx q[58], q[36];
h q[61];
h q[46];
t q[24];
cx q[79], q[24];
ccx q[38], q[62], q[20];
s q[77];
s q[26];
ccx q[56], q[48], q[28];
h q[73];
s q[47];
t q[73];
s q[48];
t q[36];
ccx q[83], q[75], q[70];
h q[82];
ccx q[36], q[66], q[21];
cx q[23], q[3];
ccx q[89], q[3], q[58];
h q[23];
cx q[66], q[25];
h q[28];
s q[82];
s q[13];
s q[46];
cx q[21], q[64];
s q[64];
ccx q[81], q[72], q[44];
cx q[52], q[56];
ccx q[6], q[16], q[31];
t q[59];
h q[13];
cx q[28], q[46];
s q[62];
s q[58];
cx q[23], q[1];
ccx q[48], q[78], q[42];
s q[30];
cx q[62], q[63];
h q[64];
h q[78];
ccx q[54], q[2], q[57];
ccx q[88], q[55], q[68];
t q[57];
ccx q[26], q[21], q[3];
t q[78];
s q[80];
cx q[51], q[9];
h q[15];
cx q[3], q[84];
cx q[64], q[42];
ccx q[89], q[47], q[13];
t q[59];
s q[57];
s q[39];
ccx q[73], q[39], q[58];
ccx q[76], q[69], q[78];
cx q[48], q[78];
h q[18];
cx q[88], q[65];
ccx q[34], q[19], q[33];
ccx q[53], q[69], q[23];
t q[79];
t q[34];
ccx q[75], q[49], q[62];
ccx q[36], q[19], q[9];
cx q[56], q[22];
s q[64];
cx q[80], q[87];
ccx q[19], q[36], q[64];
s q[28];
ccx q[87], q[11], q[85];
s q[54];
s q[22];
h q[46];
t q[10];
cx q[27], q[13];
s q[53];
cx q[67], q[54];
t q[51];
s q[31];
ccx q[29], q[6], q[54];
cx q[41], q[44];
h q[44];
h q[49];
ccx q[28], q[88], q[73];
s q[43];
h q[66];
cx q[80], q[50];
cx q[18], q[16];
t q[53];
t q[54];
t q[39];
cx q[15], q[21];
h q[70];
ccx q[63], q[89], q[18];
ccx q[38], q[84], q[13];
ccx q[23], q[5], q[45];
ccx q[55], q[59], q[70];
h q[7];
t q[75];
cx q[87], q[34];
s q[3];
s q[80];
h q[46];
t q[16];
h q[6];
ccx q[37], q[20], q[72];
h q[74];
s q[22];
h q[72];
cx q[57], q[34];
s q[23];
ccx q[57], q[26], q[12];
s q[23];
s q[40];
s q[22];
s q[76];
t q[68];
cx q[63], q[77];
cx q[36], q[28];
s q[30];
ccx q[71], q[61], q[63];
h q[12];
h q[35];
s q[84];
s q[39];
cx q[24], q[59];
t q[36];
ccx q[5], q[18], q[44];
h q[28];
h q[59];
ccx q[55], q[42], q[10];
t q[82];
t q[5];
t q[86];
ccx q[79], q[65], q[61];
cx q[4], q[16];
h q[5];
cx q[31], q[49];
cx q[53], q[65];
ccx q[12], q[56], q[86];
s q[0];
h q[54];
ccx q[72], q[69], q[17];
t q[74];
cx q[19], q[67];
h q[80];
s q[44];
s q[83];
t q[21];
s q[87];
s q[8];
ccx q[5], q[13], q[35];
ccx q[60], q[74], q[3];
h q[82];
s q[72];
s q[55];
cx q[16], q[2];
ccx q[76], q[19], q[45];
t q[30];
cx q[20], q[31];
h q[86];
h q[75];
s q[9];
t q[54];
t q[53];
s q[20];
cx q[19], q[45];
cx q[23], q[30];
ccx q[15], q[50], q[47];
ccx q[38], q[19], q[45];
h q[78];
h q[24];
h q[71];
s q[31];
s q[80];
cx q[64], q[45];
cx q[19], q[55];
ccx q[68], q[48], q[29];
cx q[11], q[77];
cx q[27], q[60];
t q[28];
h q[86];
cx q[22], q[10];
ccx q[5], q[65], q[26];
s q[71];
cx q[4], q[88];
ccx q[89], q[67], q[56];
h q[40];
h q[13];
cx q[45], q[68];
t q[43];
s q[87];
cx q[30], q[57];
cx q[39], q[63];
h q[47];
h q[10];
cx q[26], q[8];
h q[43];
ccx q[73], q[24], q[86];
s q[38];
t q[12];
s q[75];
h q[89];
t q[30];
h q[39];
t q[80];
ccx q[17], q[48], q[24];
cx q[1], q[53];
h q[32];
h q[21];
t q[45];
ccx q[21], q[14], q[10];
h q[87];
cx q[78], q[31];
s q[48];
ccx q[5], q[9], q[57];
cx q[16], q[78];
ccx q[71], q[53], q[34];
t q[30];
h q[80];
cx q[15], q[36];
ccx q[69], q[38], q[22];
cx q[40], q[45];
s q[27];
ccx q[2], q[44], q[23];
h q[61];
t q[28];
s q[0];
t q[29];
s q[33];
h q[51];
ccx q[88], q[64], q[84];
h q[15];
s q[31];
h q[25];
cx q[45], q[56];
s q[40];
h q[28];
cx q[64], q[23];
h q[61];
t q[55];
h q[21];
cx q[1], q[0];
x q[0];
z q[2];
cx q[3], q[4];
tdg q[4];
cx q[3], q[4];
tdg q[3];
sdg q[3];
sdg q[5];
tdg q[7];
sdg q[7];
tdg q[10];
sdg q[9];
y q[11];
sdg q[12];
cx q[14], q[15];
sdg q[14];
sdg q[17];
cx q[17], q[16];
tdg q[16];
sdg q[19];
tdg q[20];
cx q[19], q[20];
tdg q[19];
cx q[21], q[22];
tdg q[22];
cx q[21], q[22];
x q[21];
cx q[24], q[23];
sdg q[23];
cx q[24], q[23];
tdg q[23];
cx q[26], q[25];
tdg q[25];
cx q[26], q[25];
x q[25];
cx q[27], q[28];
x q[29];
tdg q[32];
tdg q[31];
cx q[31], q[32];
sdg q[31];
cx q[33], q[34];
tdg q[34];
cx q[33], q[34];
sdg q[34];
sdg q[33];
tdg q[35];
cx q[37], q[38];
tdg q[38];
cx q[37], q[38];
sdg q[37];
tdg q[37];
tdg q[40];
cx q[42], q[43];
sdg q[43];
cx q[42], q[43];
tdg q[43];
sdg q[42];
t q[84];
t q[46];
s q[79];
h q[77];
t q[75];
cx q[59], q[70];
cx q[65], q[47];
t q[61];
s q[50];
s q[78];
s q[66];
ccx q[78], q[60], q[52];
h q[71];
h q[51];
h q[74];
s q[73];
ccx q[76], q[78], q[88];
ccx q[47], q[65], q[52];
s q[71];
s q[88];
h q[81];
cx q[66], q[85];
t q[76];
ccx q[80], q[78], q[65];
ccx q[84], q[81], q[61];
s q[70];
s q[54];
s q[70];
ccx q[46], q[51], q[85];
ccx q[62], q[52], q[73];
cx q[79], q[88];
ccx q[49], q[63], q[73];
h q[87];
t q[85];
ccx q[55], q[51], q[89];
ccx q[52], q[84], q[81];
s q[74];
ccx q[60], q[77], q[55];
ccx q[46], q[64], q[81];
h q[89];
h q[78];
s q[55];
ccx q[84], q[59], q[47];
cx q[83], q[56];
ccx q[83], q[86], q[70];
cx q[90], q[38];
cx q[91], q[24];
cx q[92], q[10];
cx q[93], q[63];
cx q[94], q[57];
cx q[95], q[30];
cx q[96], q[28];
cx q[97], q[42];
cx q[98], q[19];
cx q[99], q[13];
cx q[100], q[18];
cx q[101], q[50];
cx q[102], q[33];
cx q[103], q[69];
cx q[104], q[71];
