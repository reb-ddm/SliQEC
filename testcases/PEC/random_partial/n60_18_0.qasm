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
ccx q[44], q[50], q[29];
s q[58];
h q[4];
s q[41];
ccx q[48], q[8], q[32];
s q[7];
s q[14];
ccx q[13], q[11], q[31];
t q[18];
h q[18];
t q[49];
cx q[43], q[57];
cx q[32], q[8];
s q[10];
h q[23];
t q[18];
ccx q[59], q[16], q[39];
cx q[1], q[18];
cx q[51], q[32];
t q[2];
t q[24];
ccx q[46], q[59], q[52];
t q[51];
ccx q[40], q[49], q[31];
ccx q[22], q[57], q[49];
t q[52];
h q[2];
h q[55];
h q[28];
ccx q[16], q[38], q[53];
t q[7];
t q[5];
s q[13];
ccx q[11], q[30], q[7];
cx q[53], q[34];
h q[45];
cx q[21], q[22];
s q[39];
t q[20];
t q[49];
cx q[22], q[28];
h q[58];
s q[21];
ccx q[28], q[26], q[35];
t q[14];
h q[58];
cx q[32], q[10];
s q[4];
s q[28];
ccx q[8], q[58], q[52];
cx q[32], q[54];
h q[34];
s q[25];
ccx q[3], q[16], q[55];
s q[9];
t q[25];
cx q[30], q[23];
ccx q[44], q[1], q[16];
ccx q[40], q[3], q[16];
cx q[17], q[29];
s q[26];
cx q[24], q[18];
t q[30];
ccx q[54], q[16], q[57];
s q[22];
t q[23];
t q[25];
ccx q[12], q[42], q[25];
ccx q[37], q[48], q[27];
t q[39];
cx q[41], q[3];
ccx q[41], q[26], q[5];
h q[45];
h q[51];
s q[58];
s q[23];
cx q[47], q[0];
h q[11];
h q[5];
cx q[48], q[15];
t q[21];
t q[23];
h q[50];
t q[13];
cx q[32], q[40];
ccx q[9], q[58], q[53];
h q[41];
h q[40];
s q[22];
t q[47];
cx q[36], q[12];
cx q[33], q[18];
ccx q[16], q[24], q[25];
s q[41];
t q[34];
t q[46];
h q[35];
h q[4];
s q[31];
cx q[52], q[34];
s q[20];
s q[24];
cx q[35], q[20];
ccx q[30], q[53], q[40];
h q[9];
h q[3];
h q[54];
ccx q[46], q[26], q[4];
t q[39];
cx q[59], q[1];
cx q[54], q[18];
s q[3];
cx q[52], q[46];
ccx q[27], q[0], q[6];
h q[21];
s q[26];
h q[6];
cx q[30], q[48];
s q[3];
t q[15];
cx q[43], q[49];
t q[9];
cx q[9], q[39];
cx q[33], q[29];
cx q[28], q[6];
cx q[34], q[12];
s q[30];
h q[36];
cx q[25], q[1];
h q[55];
h q[0];
t q[47];
ccx q[29], q[24], q[50];
h q[39];
cx q[16], q[45];
h q[49];
cx q[42], q[53];
s q[54];
cx q[18], q[44];
cx q[46], q[3];
s q[5];
cx q[13], q[49];
t q[37];
ccx q[50], q[23], q[22];
t q[3];
h q[25];
t q[35];
cx q[28], q[8];
h q[13];
ccx q[33], q[24], q[17];
s q[53];
s q[48];
ccx q[49], q[17], q[12];
t q[47];
cx q[4], q[58];
s q[7];
s q[28];
h q[14];
cx q[10], q[24];
t q[58];
ccx q[4], q[50], q[59];
s q[19];
h q[43];
h q[12];
h q[50];
ccx q[5], q[12], q[33];
ccx q[38], q[29], q[28];
t q[36];
h q[27];
h q[32];
h q[18];
cx q[39], q[31];
ccx q[19], q[47], q[28];
t q[33];
s q[55];
h q[32];
s q[19];
s q[42];
ccx q[48], q[44], q[47];
s q[0];
tdg q[1];
sdg q[0];
cx q[0], q[1];
x q[2];
cx q[4], q[3];
tdg q[3];
cx q[5], q[6];
sdg q[6];
cx q[5], q[6];
tdg q[5];
x q[7];
cx q[9], q[10];
tdg q[9];
x q[11];
x q[12];
sdg q[13];
sdg q[15];
cx q[18], q[17];
sdg q[17];
cx q[19], q[20];
sdg q[20];
cx q[19], q[20];
tdg q[20];
sdg q[19];
x q[22];
sdg q[23];
cx q[23], q[24];
sdg q[26];
cx q[25], q[26];
sdg q[25];
cx q[27], q[28];
sdg q[28];
cx q[27], q[28];
tdg q[27];
h q[42];
ccx q[30], q[38], q[56];
t q[39];
ccx q[35], q[38], q[54];
cx q[47], q[57];
h q[34];
t q[50];
t q[49];
s q[34];
h q[46];
h q[43];
t q[55];
t q[54];
t q[34];
t q[34];
h q[36];
t q[49];
s q[44];
h q[59];
ccx q[40], q[56], q[55];
h q[39];
cx q[38], q[31];
ccx q[44], q[56], q[48];
s q[57];
ccx q[33], q[52], q[51];
ccx q[40], q[47], q[37];
h q[57];
t q[33];
cx q[32], q[30];
cx q[40], q[44];
