OPENQASM 2.0;
include "qelib1.inc";
qreg q[40];
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
s q[36];
s q[13];
s q[12];
s q[13];
h q[15];
s q[29];
h q[34];
t q[29];
t q[9];
s q[10];
t q[16];
ccx q[18], q[34], q[19];
s q[9];
s q[1];
h q[28];
s q[31];
s q[7];
t q[23];
cx q[12], q[5];
t q[12];
s q[35];
cx q[4], q[29];
ccx q[29], q[18], q[2];
cx q[22], q[1];
t q[14];
s q[5];
s q[9];
cx q[28], q[7];
h q[38];
cx q[19], q[37];
ccx q[16], q[37], q[39];
ccx q[31], q[27], q[36];
s q[14];
cx q[3], q[6];
ccx q[20], q[10], q[7];
cx q[4], q[30];
s q[25];
s q[3];
ccx q[14], q[32], q[11];
h q[39];
cx q[11], q[10];
t q[17];
t q[27];
s q[18];
s q[33];
ccx q[38], q[27], q[5];
s q[11];
s q[33];
t q[18];
h q[39];
cx q[13], q[18];
t q[38];
cx q[11], q[39];
s q[32];
h q[5];
ccx q[39], q[34], q[27];
h q[17];
t q[2];
h q[30];
cx q[25], q[14];
cx q[39], q[23];
s q[26];
t q[23];
t q[4];
t q[10];
h q[23];
ccx q[33], q[14], q[24];
h q[6];
s q[9];
cx q[10], q[13];
h q[31];
s q[12];
cx q[11], q[0];
cx q[32], q[3];
ccx q[1], q[6], q[30];
s q[36];
ccx q[29], q[1], q[11];
t q[30];
ccx q[17], q[33], q[13];
t q[39];
cx q[32], q[4];
h q[36];
ccx q[5], q[27], q[23];
ccx q[33], q[14], q[3];
h q[28];
cx q[30], q[32];
t q[39];
t q[11];
s q[2];
s q[10];
t q[2];
t q[36];
cx q[27], q[10];
h q[10];
cx q[4], q[39];
t q[34];
cx q[36], q[21];
ccx q[10], q[9], q[0];
cx q[11], q[15];
h q[20];
s q[31];
ccx q[36], q[32], q[18];
ccx q[12], q[4], q[13];
t q[13];
s q[34];
h q[16];
h q[30];
cx q[38], q[37];
t q[32];
ccx q[12], q[1], q[32];
t q[17];
h q[19];
cx q[0], q[33];
cx q[5], q[10];
s q[14];
s q[20];
ccx q[24], q[14], q[25];
cx q[20], q[22];
s q[18];
t q[23];
