OPENQASM 2.0;
include "qelib1.inc";
qreg q[117];
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
t q[37];
cx q[93], q[3];
cx q[72], q[13];
h q[21];
cx q[48], q[7];
t q[15];
s q[86];
ccx q[26], q[75], q[70];
h q[86];
s q[68];
t q[24];
t q[2];
ccx q[19], q[20], q[16];
h q[87];
ccx q[64], q[51], q[68];
t q[98];
t q[28];
ccx q[56], q[65], q[71];
ccx q[13], q[21], q[92];
ccx q[20], q[24], q[47];
t q[87];
cx q[11], q[80];
s q[92];
h q[33];
h q[17];
t q[98];
h q[68];
ccx q[32], q[65], q[24];
h q[87];
ccx q[75], q[87], q[30];
t q[53];
t q[24];
h q[55];
h q[8];
cx q[25], q[92];
h q[93];
h q[45];
s q[65];
ccx q[60], q[62], q[13];
cx q[97], q[76];
t q[19];
t q[83];
ccx q[14], q[96], q[45];
h q[75];
t q[72];
t q[3];
t q[34];
h q[15];
t q[55];
h q[88];
t q[77];
h q[39];
s q[69];
ccx q[96], q[8], q[87];
s q[10];
ccx q[5], q[8], q[66];
s q[6];
t q[34];
t q[8];
t q[80];
cx q[73], q[12];
ccx q[55], q[58], q[78];
s q[41];
h q[86];
t q[19];
t q[60];
t q[46];
h q[6];
ccx q[7], q[98], q[14];
t q[18];
t q[93];
s q[34];
h q[91];
s q[21];
t q[30];
ccx q[64], q[96], q[79];
h q[74];
t q[47];
t q[34];
ccx q[49], q[43], q[21];
t q[68];
h q[48];
h q[95];
cx q[92], q[3];
cx q[87], q[41];
ccx q[37], q[75], q[28];
h q[24];
h q[39];
s q[80];
h q[5];
ccx q[16], q[57], q[63];
ccx q[53], q[75], q[71];
h q[38];
cx q[55], q[9];
h q[69];
h q[87];
t q[85];
t q[58];
cx q[89], q[29];
s q[37];
ccx q[42], q[66], q[37];
s q[31];
h q[65];
cx q[37], q[15];
cx q[13], q[3];
h q[30];
s q[98];
ccx q[45], q[20], q[62];
ccx q[65], q[61], q[54];
t q[10];
t q[54];
cx q[47], q[17];
ccx q[79], q[12], q[13];
s q[38];
s q[42];
s q[86];
s q[46];
ccx q[22], q[39], q[85];
ccx q[14], q[10], q[23];
ccx q[12], q[6], q[88];
cx q[34], q[94];
ccx q[30], q[26], q[48];
t q[80];
s q[2];
ccx q[44], q[7], q[38];
h q[52];
t q[31];
s q[75];
h q[33];
t q[80];
h q[3];
ccx q[86], q[96], q[29];
t q[15];
cx q[35], q[56];
h q[19];
ccx q[98], q[71], q[80];
ccx q[24], q[52], q[19];
s q[48];
s q[65];
t q[45];
s q[20];
ccx q[61], q[56], q[27];
h q[75];
t q[72];
ccx q[77], q[71], q[45];
ccx q[85], q[38], q[63];
ccx q[95], q[15], q[46];
s q[44];
t q[76];
h q[36];
s q[63];
s q[87];
h q[21];
h q[74];
h q[96];
t q[49];
t q[22];
ccx q[41], q[1], q[31];
h q[70];
cx q[76], q[28];
h q[55];
s q[9];
ccx q[9], q[47], q[25];
ccx q[96], q[26], q[21];
t q[55];
t q[71];
h q[44];
s q[78];
s q[58];
s q[66];
h q[31];
h q[0];
t q[75];
t q[69];
t q[60];
s q[74];
t q[79];
h q[32];
t q[64];
cx q[38], q[53];
t q[96];
s q[60];
h q[13];
t q[57];
t q[36];
s q[55];
cx q[43], q[5];
ccx q[90], q[0], q[48];
cx q[46], q[54];
h q[83];
ccx q[81], q[7], q[39];
t q[62];
s q[35];
ccx q[51], q[72], q[87];
s q[45];
s q[10];
cx q[46], q[66];
cx q[19], q[91];
h q[28];
s q[74];
s q[37];
cx q[92], q[65];
t q[61];
ccx q[65], q[10], q[89];
t q[39];
cx q[63], q[97];
h q[8];
h q[72];
s q[48];
s q[17];
s q[30];
t q[73];
ccx q[16], q[44], q[13];
t q[58];
t q[30];
t q[40];
ccx q[81], q[65], q[61];
ccx q[57], q[56], q[42];
cx q[78], q[50];
s q[68];
h q[46];
s q[73];
cx q[94], q[33];
h q[5];
cx q[52], q[12];
cx q[20], q[43];
h q[32];
t q[58];
cx q[36], q[67];
h q[82];
s q[69];
t q[22];
s q[66];
ccx q[93], q[51], q[31];
t q[62];
t q[43];
cx q[10], q[95];
t q[18];
s q[34];
h q[71];
cx q[64], q[38];
s q[86];
ccx q[32], q[74], q[22];
s q[34];
s q[89];
s q[73];
s q[48];
t q[47];
cx q[18], q[62];
h q[35];
h q[8];
t q[40];
ccx q[29], q[88], q[93];
h q[81];
s q[49];
cx q[49], q[3];
cx q[13], q[74];
ccx q[11], q[85], q[2];
s q[71];
t q[5];
t q[63];
ccx q[99], q[35], q[27];
t q[80];
ccx q[36], q[90], q[83];
t q[31];
t q[37];
t q[74];
h q[81];
ccx q[85], q[35], q[54];
h q[6];
h q[91];
t q[70];
s q[62];
cx q[91], q[62];
cx q[36], q[25];
cx q[55], q[40];
ccx q[91], q[12], q[54];
h q[18];
cx q[45], q[24];
h q[13];
cx q[56], q[34];
cx q[79], q[93];
cx q[98], q[77];
s q[52];
s q[91];
h q[88];
s q[51];
s q[39];
t q[65];
t q[31];
t q[82];
ccx q[74], q[93], q[85];
s q[6];
cx q[0], q[80];
t q[16];
ccx q[49], q[63], q[40];
h q[96];
t q[61];
s q[92];
t q[46];
sdg q[0];
y q[2];
cx q[4], q[3];
tdg q[3];
x q[3];
cx q[4], q[3];
x q[3];
sdg q[5];
tdg q[5];
sdg q[7];
tdg q[8];
cx q[8], q[7];
tdg q[7];
y q[9];
tdg q[11];
cx q[10], q[11];
sdg q[10];
x q[12];
y q[13];
tdg q[14];
cx q[15], q[14];
sdg q[14];
sdg q[16];
cx q[17], q[16];
tdg q[16];
tdg q[18];
cx q[19], q[18];
tdg q[18];
cx q[19], q[18];
sdg q[18];
tdg q[20];
z q[22];
z q[23];
x q[24];
cx q[24], q[25];
sdg q[25];
cx q[24], q[25];
x q[24];
cx q[26], q[27];
tdg q[27];
sdg q[26];
cx q[26], q[27];
sdg q[26];
sdg q[28];
cx q[28], q[29];
tdg q[29];
tdg q[28];
cx q[28], q[29];
cx q[30], q[31];
sdg q[31];
cx q[30], q[31];
x q[30];
z q[32];
x q[33];
x q[34];
z q[35];
y q[36];
tdg q[37];
sdg q[39];
cx q[39], q[40];
sdg q[40];
cx q[39], q[40];
sdg q[39];
cx q[41], q[42];
sdg q[41];
cx q[44], q[43];
x q[43];
tdg q[45];
tdg q[47];
x q[47];
cx q[48], q[47];
x q[47];
y q[49];
cx q[75], q[84];
s q[92];
h q[82];
t q[80];
t q[71];
t q[98];
t q[73];
ccx q[59], q[77], q[56];
h q[87];
t q[74];
ccx q[72], q[87], q[91];
h q[61];
h q[84];
cx q[87], q[52];
cx q[51], q[73];
ccx q[74], q[54], q[91];
s q[85];
t q[72];
cx q[76], q[81];
cx q[93], q[67];
ccx q[75], q[76], q[58];
s q[92];
t q[76];
cx q[65], q[56];
cx q[95], q[65];
ccx q[79], q[57], q[67];
h q[93];
ccx q[96], q[53], q[72];
t q[54];
t q[88];
s q[92];
ccx q[76], q[92], q[94];
cx q[59], q[78];
ccx q[92], q[63], q[71];
h q[53];
cx q[72], q[89];
s q[86];
h q[88];
s q[57];
t q[91];
h q[79];
h q[95];
cx q[53], q[67];
s q[86];
h q[94];
ccx q[81], q[59], q[86];
s q[77];
ccx q[54], q[78], q[70];
s q[84];
h q[89];
cx q[100], q[22];
cx q[101], q[47];
cx q[102], q[92];
cx q[103], q[6];
cx q[104], q[33];
cx q[105], q[74];
cx q[106], q[69];
cx q[107], q[61];
cx q[108], q[98];
cx q[109], q[30];
cx q[110], q[48];
cx q[111], q[81];
cx q[112], q[55];
cx q[113], q[11];
cx q[114], q[62];
cx q[115], q[60];
cx q[116], q[66];
