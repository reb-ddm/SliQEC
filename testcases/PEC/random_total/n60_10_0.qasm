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
ccx q[22], q[40], q[55];
ccx q[45], q[4], q[10];
s q[53];
s q[37];
s q[17];
h q[25];
ccx q[0], q[25], q[14];
h q[42];
cx q[23], q[5];
ccx q[43], q[39], q[54];
cx q[25], q[28];
t q[0];
t q[13];
s q[31];
t q[48];
t q[53];
h q[19];
t q[55];
cx q[48], q[58];
s q[56];
h q[46];
h q[5];
s q[57];
cx q[11], q[1];
ccx q[36], q[25], q[11];
s q[25];
ccx q[32], q[14], q[39];
s q[42];
ccx q[19], q[33], q[43];
h q[17];
t q[23];
s q[5];
h q[44];
cx q[23], q[29];
h q[58];
h q[34];
cx q[34], q[21];
cx q[9], q[55];
s q[56];
s q[26];
cx q[58], q[39];
h q[9];
s q[3];
s q[2];
s q[20];
ccx q[32], q[59], q[15];
t q[14];
cx q[42], q[8];
ccx q[9], q[8], q[2];
cx q[17], q[29];
t q[50];
h q[8];
cx q[48], q[30];
ccx q[56], q[6], q[10];
cx q[12], q[27];
ccx q[37], q[41], q[56];
s q[7];
s q[5];
t q[28];
t q[2];
cx q[44], q[46];
cx q[17], q[31];
cx q[59], q[3];
s q[20];
ccx q[55], q[53], q[7];
t q[46];
t q[20];
cx q[46], q[54];
ccx q[45], q[9], q[52];
s q[22];
t q[32];
t q[29];
t q[1];
h q[47];
ccx q[31], q[51], q[23];
h q[47];
cx q[13], q[55];
t q[19];
s q[35];
s q[45];
t q[42];
h q[44];
cx q[20], q[52];
h q[19];
ccx q[12], q[42], q[10];
cx q[3], q[22];
cx q[51], q[26];
s q[39];
h q[53];
t q[29];
s q[33];
t q[38];
t q[29];
h q[48];
s q[23];
t q[22];
cx q[59], q[18];
t q[41];
s q[53];
t q[55];
ccx q[3], q[44], q[27];
ccx q[51], q[59], q[1];
h q[38];
h q[47];
h q[40];
cx q[39], q[12];
h q[8];
t q[42];
t q[24];
cx q[31], q[14];
ccx q[40], q[36], q[6];
ccx q[38], q[21], q[48];
s q[0];
ccx q[42], q[34], q[7];
s q[10];
cx q[8], q[45];
cx q[54], q[45];
h q[58];
h q[15];
cx q[25], q[20];
h q[57];
s q[22];
ccx q[2], q[57], q[9];
t q[36];
s q[0];
h q[26];
t q[57];
ccx q[3], q[16], q[10];
h q[46];
t q[32];
t q[14];
ccx q[50], q[34], q[37];
t q[6];
t q[18];
h q[12];
cx q[1], q[39];
s q[59];
t q[53];
cx q[37], q[22];
h q[5];
h q[27];
ccx q[3], q[48], q[49];
cx q[37], q[42];
s q[55];
s q[2];
h q[8];
ccx q[26], q[39], q[8];
s q[9];
t q[40];
h q[27];
t q[6];
cx q[20], q[25];
t q[26];
cx q[43], q[23];
s q[52];
s q[15];
cx q[54], q[4];
t q[27];
h q[39];
h q[20];
cx q[37], q[20];
h q[6];
cx q[44], q[38];
s q[11];
t q[35];
h q[24];
t q[30];
s q[49];
s q[53];
t q[30];
h q[8];
h q[58];
ccx q[0], q[48], q[15];
s q[18];
h q[17];
ccx q[36], q[54], q[20];
ccx q[55], q[27], q[33];
cx q[39], q[29];
s q[38];
s q[7];