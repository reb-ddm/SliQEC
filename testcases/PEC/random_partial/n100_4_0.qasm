OPENQASM 2.0;
include "qelib1.inc";
qreg q[100];
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
h q[90];
h q[91];
h q[92];
h q[93];
h q[94];
h q[95];
h q[96];
h q[97];
h q[98];
h q[99];
h q[56];
h q[62];
s q[9];
cx q[65], q[58];
t q[12];
cx q[12], q[35];
cx q[57], q[19];
s q[26];
ccx q[42], q[81], q[15];
cx q[78], q[21];
h q[59];
s q[13];
ccx q[42], q[93], q[69];
t q[46];
cx q[22], q[23];
cx q[10], q[63];
h q[38];
h q[50];
ccx q[84], q[28], q[88];
t q[64];
s q[68];
cx q[85], q[63];
cx q[70], q[75];
cx q[21], q[89];
h q[33];
cx q[3], q[16];
s q[14];
s q[97];
t q[45];
t q[16];
t q[78];
ccx q[38], q[4], q[25];
ccx q[92], q[14], q[69];
ccx q[24], q[30], q[16];
h q[86];
h q[10];
t q[7];
h q[4];
ccx q[6], q[18], q[66];
h q[29];
s q[34];
ccx q[1], q[98], q[21];
cx q[68], q[92];
t q[20];
cx q[24], q[68];
ccx q[39], q[73], q[48];
h q[89];
ccx q[12], q[43], q[52];
ccx q[18], q[40], q[90];
cx q[62], q[45];
ccx q[45], q[1], q[21];
t q[13];
cx q[30], q[31];
ccx q[17], q[77], q[60];
s q[83];
t q[76];
h q[82];
cx q[46], q[50];
cx q[67], q[54];
cx q[7], q[63];
cx q[65], q[33];
cx q[46], q[87];
t q[2];
h q[4];
t q[81];
h q[6];
ccx q[84], q[81], q[73];
cx q[13], q[99];
t q[97];
cx q[7], q[29];
s q[37];
cx q[12], q[74];
s q[21];
ccx q[50], q[19], q[56];
h q[77];
cx q[63], q[1];
ccx q[32], q[81], q[89];
cx q[6], q[61];
s q[63];
cx q[84], q[80];
t q[48];
t q[63];
h q[43];
ccx q[40], q[54], q[44];
s q[73];
s q[88];
t q[31];
cx q[23], q[72];
cx q[97], q[85];
h q[54];
ccx q[26], q[62], q[8];
cx q[60], q[16];
t q[4];
cx q[32], q[68];
t q[46];
h q[68];
t q[53];
ccx q[16], q[74], q[0];
t q[30];
h q[94];
cx q[7], q[44];
s q[33];
s q[22];
h q[77];
h q[80];
h q[26];
cx q[11], q[38];
h q[33];
cx q[29], q[9];
cx q[39], q[97];
ccx q[79], q[48], q[37];
h q[50];
ccx q[86], q[91], q[32];
t q[12];
ccx q[35], q[12], q[30];
cx q[72], q[44];
ccx q[0], q[17], q[88];
t q[95];
t q[28];
h q[96];
h q[96];
t q[91];
ccx q[13], q[14], q[46];
t q[95];
cx q[72], q[42];
cx q[88], q[85];
h q[4];
h q[42];
h q[44];
h q[32];
s q[87];
t q[72];
t q[40];
ccx q[93], q[24], q[89];
h q[67];
t q[78];
h q[60];
h q[62];
cx q[46], q[29];
h q[48];
h q[34];
s q[78];
cx q[94], q[44];
t q[98];
ccx q[53], q[15], q[38];
h q[84];
cx q[93], q[98];
h q[89];
ccx q[43], q[31], q[6];
t q[35];
ccx q[43], q[96], q[46];
ccx q[52], q[30], q[2];
h q[49];
t q[90];
ccx q[1], q[73], q[91];
s q[2];
cx q[17], q[6];
ccx q[96], q[77], q[31];
ccx q[27], q[79], q[32];
t q[70];
s q[54];
h q[51];
t q[45];
s q[65];
t q[72];
t q[24];
s q[14];
h q[82];
s q[18];
ccx q[1], q[76], q[12];
s q[15];
s q[96];
h q[61];
cx q[14], q[64];
s q[96];
h q[93];
cx q[71], q[32];
s q[62];
t q[28];
h q[52];
ccx q[27], q[14], q[97];
ccx q[55], q[28], q[69];
t q[0];
cx q[45], q[52];
s q[54];
ccx q[78], q[46], q[63];
cx q[91], q[42];
s q[75];
cx q[88], q[21];
cx q[77], q[93];
cx q[29], q[90];
cx q[34], q[61];
t q[50];
h q[9];
s q[52];
cx q[55], q[63];
t q[30];
ccx q[22], q[2], q[81];
s q[67];
t q[83];
cx q[11], q[72];
cx q[2], q[48];
t q[35];
h q[72];
s q[49];
ccx q[26], q[31], q[67];
t q[28];
t q[41];
s q[41];
ccx q[44], q[96], q[31];
t q[88];
ccx q[97], q[38], q[80];
s q[15];
s q[45];
h q[62];
t q[17];
s q[99];
ccx q[89], q[57], q[20];
ccx q[74], q[84], q[93];
cx q[51], q[57];
t q[49];
ccx q[95], q[35], q[63];
s q[53];
h q[78];
cx q[47], q[46];
cx q[6], q[3];
t q[57];
cx q[84], q[62];
cx q[14], q[75];
h q[0];
h q[1];
ccx q[93], q[98], q[60];
t q[62];
ccx q[90], q[97], q[16];
s q[66];
h q[57];
s q[91];
h q[14];
ccx q[71], q[28], q[97];
ccx q[18], q[8], q[30];
s q[50];
t q[99];
h q[4];
ccx q[0], q[44], q[68];
h q[77];
ccx q[23], q[36], q[54];
s q[1];
cx q[48], q[9];
s q[26];
cx q[81], q[19];
h q[98];
ccx q[94], q[56], q[60];
ccx q[81], q[93], q[90];
ccx q[73], q[47], q[86];
t q[89];
s q[45];
t q[9];
cx q[14], q[41];
ccx q[61], q[73], q[29];
t q[83];
cx q[62], q[39];
ccx q[45], q[88], q[25];
t q[81];
h q[23];
cx q[26], q[48];
s q[20];
h q[93];
t q[1];
h q[45];
s q[41];
s q[40];
h q[7];
cx q[32], q[92];
h q[48];
ccx q[69], q[47], q[61];
cx q[21], q[72];
cx q[82], q[12];
h q[8];
s q[21];
ccx q[75], q[68], q[67];
s q[32];
ccx q[99], q[7], q[56];
cx q[85], q[95];
ccx q[69], q[76], q[72];
t q[57];
cx q[78], q[9];
cx q[86], q[74];
ccx q[85], q[80], q[5];
s q[62];
t q[92];
cx q[53], q[24];
s q[64];
h q[19];
ccx q[46], q[56], q[47];
cx q[68], q[66];
t q[69];
h q[61];
ccx q[7], q[92], q[63];
t q[83];
h q[12];
x q[0];
cx q[2], q[3];
tdg q[2];
cx q[5], q[4];
sdg q[4];
sdg q[5];
cx q[5], q[4];
sdg q[4];
cx q[6], q[7];
sdg q[6];
tdg q[7];
cx q[6], q[7];
sdg q[6];
z q[8];
y q[9];
y q[10];
z q[11];
sdg q[12];
cx q[14], q[15];
sdg q[14];
cx q[16], q[17];
sdg q[17];
cx q[16], q[17];
sdg q[17];
tdg q[16];
y q[18];
x q[19];
tdg q[21];
cx q[22], q[21];
tdg q[21];
cx q[22], q[21];
sdg q[21];
sdg q[24];
tdg q[24];
cx q[23], q[24];
sdg q[23];
sdg q[26];
sdg q[25];
cx q[25], q[26];
tdg q[27];
z q[29];
tdg q[31];
cx q[30], q[31];
sdg q[30];
tdg q[30];
tdg q[33];
cx q[33], q[34];
cx q[35], q[36];
tdg q[35];
cx q[38], q[37];
sdg q[37];
cx q[38], q[37];
tdg q[37];
sdg q[39];
sdg q[40];
cx q[40], q[39];
sdg q[39];
sdg q[41];
tdg q[44];
cx q[44], q[45];
tdg q[45];
sdg q[44];
cx q[44], q[45];
cx q[47], q[46];
tdg q[47];
tdg q[46];
cx q[47], q[46];
sdg q[46];
cx q[48], q[49];
tdg q[49];
cx q[48], q[49];
tdg q[48];
s q[56];
cx q[89], q[66];
h q[96];
ccx q[59], q[71], q[94];
ccx q[61], q[93], q[63];
cx q[64], q[79];
s q[83];
s q[53];
ccx q[56], q[80], q[50];
cx q[62], q[53];
cx q[74], q[53];
cx q[59], q[92];
t q[54];
ccx q[69], q[88], q[53];
ccx q[77], q[56], q[74];
cx q[54], q[75];
s q[80];
t q[64];
cx q[65], q[58];
t q[89];
cx q[71], q[54];
s q[95];
cx q[90], q[96];
ccx q[52], q[66], q[72];
h q[59];
cx q[52], q[80];
h q[92];
cx q[58], q[65];
ccx q[69], q[55], q[72];
ccx q[81], q[55], q[70];
h q[98];
cx q[78], q[95];
h q[80];
ccx q[93], q[81], q[76];
ccx q[82], q[76], q[58];
t q[74];
cx q[81], q[96];
t q[54];
ccx q[61], q[87], q[98];
cx q[74], q[84];
t q[58];
cx q[85], q[70];
cx q[56], q[93];
h q[57];
ccx q[95], q[56], q[96];
s q[69];
t q[70];
s q[79];
h q[99];
cx q[64], q[68];
