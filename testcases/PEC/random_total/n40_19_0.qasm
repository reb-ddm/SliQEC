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
h q[35];
cx q[37], q[25];
cx q[11], q[20];
ccx q[36], q[13], q[28];
cx q[33], q[13];
h q[0];
ccx q[1], q[15], q[24];
h q[8];
h q[11];
t q[23];
s q[3];
cx q[21], q[29];
ccx q[39], q[15], q[32];
t q[20];
t q[23];
s q[37];
h q[1];
cx q[3], q[21];
h q[13];
s q[11];
ccx q[36], q[6], q[35];
h q[22];
t q[8];
ccx q[17], q[16], q[2];
t q[24];
cx q[20], q[9];
s q[32];
ccx q[9], q[27], q[12];
s q[5];
ccx q[37], q[36], q[8];
cx q[32], q[17];
h q[11];
ccx q[1], q[9], q[19];
cx q[15], q[11];
ccx q[38], q[23], q[33];
h q[15];
ccx q[27], q[34], q[31];
h q[36];
h q[6];
s q[3];
h q[23];
s q[3];
cx q[3], q[35];
s q[16];
s q[21];
t q[14];
cx q[0], q[36];
ccx q[16], q[28], q[21];
cx q[26], q[20];
t q[2];
h q[1];
t q[30];
cx q[11], q[15];
h q[4];
s q[16];
cx q[0], q[7];
ccx q[31], q[17], q[35];
ccx q[39], q[19], q[14];
h q[10];
s q[39];
h q[10];
h q[38];
h q[36];
ccx q[18], q[33], q[29];
cx q[0], q[8];
s q[36];
s q[30];
cx q[39], q[36];
cx q[16], q[15];
ccx q[33], q[1], q[25];
cx q[9], q[21];
h q[15];
h q[13];
t q[27];
s q[11];
h q[5];
cx q[3], q[21];
h q[31];
h q[25];
ccx q[6], q[14], q[22];
ccx q[33], q[18], q[26];
t q[16];
s q[2];
cx q[20], q[7];
t q[36];
t q[25];
s q[5];
s q[22];
s q[35];
ccx q[2], q[20], q[4];
cx q[33], q[34];
h q[25];
s q[35];
h q[31];
s q[23];
t q[0];
cx q[24], q[2];
h q[24];
ccx q[3], q[15], q[17];
ccx q[16], q[29], q[18];
h q[13];
cx q[7], q[29];
t q[10];
cx q[17], q[4];
t q[30];
ccx q[15], q[35], q[5];
t q[12];
h q[30];
h q[2];
t q[12];
h q[24];
t q[31];
h q[35];
t q[11];
t q[3];
h q[1];
cx q[13], q[38];
t q[26];
s q[1];
s q[15];
