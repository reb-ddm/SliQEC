OPENQASM 2.0;
include "qelib1.inc";
qreg q[60];
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
h q[45];
ccx q[51], q[4], q[44];
t q[59];
h q[41];
ccx q[59], q[10], q[57];
ccx q[46], q[28], q[24];
cx q[42], q[53];
h q[27];
ccx q[46], q[36], q[3];
cx q[13], q[54];
t q[56];
ccx q[53], q[59], q[4];
ccx q[3], q[18], q[1];
h q[56];
cx q[37], q[4];
cx q[20], q[48];
h q[27];
s q[11];
ccx q[6], q[44], q[48];
cx q[53], q[44];
s q[50];
ccx q[14], q[5], q[46];
s q[36];
s q[47];
ccx q[43], q[36], q[59];
ccx q[26], q[25], q[59];
s q[17];
cx q[56], q[22];
t q[41];
t q[52];
cx q[50], q[9];
t q[24];
t q[25];
t q[21];
cx q[14], q[45];
h q[27];
h q[22];
h q[48];
cx q[57], q[52];
s q[21];
ccx q[10], q[12], q[39];
t q[31];
h q[59];
cx q[18], q[21];
t q[51];
h q[25];
cx q[39], q[33];
h q[15];
t q[54];
t q[13];
t q[7];
t q[43];
t q[47];
h q[49];
cx q[13], q[48];
h q[52];
t q[49];
cx q[50], q[17];
ccx q[52], q[31], q[29];
t q[12];
h q[33];
cx q[51], q[22];
t q[6];
t q[31];
cx q[52], q[0];
s q[34];
s q[45];
cx q[39], q[57];
ccx q[26], q[51], q[35];
cx q[53], q[2];
ccx q[38], q[1], q[9];
h q[32];
h q[8];
cx q[5], q[6];
cx q[38], q[17];
t q[13];
ccx q[49], q[43], q[1];
t q[36];
ccx q[41], q[22], q[8];
h q[31];
ccx q[35], q[46], q[41];
h q[21];
cx q[8], q[24];
h q[15];
h q[4];
t q[39];
h q[28];
ccx q[28], q[54], q[23];
s q[1];
s q[7];
cx q[31], q[34];
cx q[4], q[52];
s q[16];
cx q[27], q[13];
cx q[21], q[11];
s q[40];
cx q[46], q[53];
ccx q[29], q[35], q[52];
h q[19];
s q[39];
cx q[52], q[31];
cx q[42], q[48];
t q[56];
t q[53];
t q[45];
ccx q[26], q[13], q[16];
cx q[41], q[32];
h q[27];
t q[55];
t q[25];
t q[56];
cx q[35], q[40];
cx q[16], q[19];
h q[14];
cx q[13], q[0];
cx q[8], q[2];
h q[55];
t q[41];
ccx q[31], q[39], q[21];
t q[35];
h q[19];
ccx q[20], q[5], q[34];
t q[1];
ccx q[11], q[55], q[3];
h q[36];
s q[43];
h q[11];
ccx q[56], q[4], q[13];
ccx q[31], q[21], q[59];
cx q[22], q[47];
h q[32];
cx q[53], q[16];
t q[39];
h q[16];
s q[49];
cx q[12], q[32];
t q[53];
h q[5];
cx q[26], q[56];
t q[58];
h q[52];
s q[18];
cx q[49], q[2];
ccx q[53], q[17], q[34];
t q[58];
h q[50];
ccx q[2], q[16], q[52];
cx q[29], q[25];
s q[6];
cx q[58], q[54];
cx q[20], q[23];
t q[1];
ccx q[12], q[36], q[13];
h q[12];
cx q[17], q[30];
h q[10];
h q[19];
cx q[56], q[16];
t q[48];
ccx q[39], q[14], q[59];
h q[48];
cx q[50], q[7];
cx q[28], q[11];
s q[14];
h q[54];
s q[38];
h q[29];
h q[12];
t q[1];
s q[4];
t q[42];
h q[3];
cx q[6], q[27];
t q[10];
cx q[38], q[12];
cx q[12], q[29];
h q[28];
t q[13];
t q[48];
cx q[52], q[4];
tdg q[0];
tdg q[3];
sdg q[2];
cx q[2], q[3];
sdg q[4];
sdg q[6];
cx q[6], q[7];
tdg q[8];
cx q[9], q[8];
tdg q[8];
sdg q[10];
sdg q[11];
cx q[10], q[11];
tdg q[10];
x q[12];
cx q[14], q[13];
sdg q[13];
tdg q[16];
sdg q[15];
cx q[15], q[16];
y q[19];
x q[20];
cx q[22], q[21];
tdg q[22];
tdg q[21];
cx q[24], q[25];
tdg q[24];
sdg q[27];
cx q[27], q[28];
z q[29];
cx q[47], q[40];
ccx q[53], q[31], q[46];
s q[38];
h q[37];
ccx q[34], q[59], q[36];
cx q[42], q[48];
t q[55];
t q[52];
s q[36];
cx q[51], q[40];
s q[33];
t q[37];
t q[33];
ccx q[46], q[54], q[42];
ccx q[36], q[50], q[48];
h q[49];
t q[32];
ccx q[56], q[43], q[34];
cx q[57], q[38];
t q[39];
ccx q[57], q[41], q[56];
cx q[37], q[51];
h q[51];
s q[54];
t q[39];
cx q[54], q[38];
t q[39];
h q[33];
t q[55];
h q[44];
