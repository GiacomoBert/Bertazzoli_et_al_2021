% MSWCMPTP  1 ���������M���̈��k�̃X���b�V���z�[���h�Ɛ��\
%
%   [THR_VAL,L2_Perf,N0_Perf] = MSWCMPTP(DEC,METH) �܂��� 
%   [THR_VAL,L2_Perf,N0_Perf] = MSWCMPTP(DEC,METH,PARAM) �́A
%   METH ��@�ƁA�K�v�ȏꍇ PARAM �p�����[�^���g���āA���k��ɓ�����
%   �x�N�g�� THR_VAL, L2_Perf, N0_Perf ���v�Z���܂� (METH �� PARAM ��
%   ���ڂ������ɂ��Ă� MSWCMP ���Q��)�B
%   i �Ԗڂ̐M���ɑ΂��āA
%     - THR_VAL(i) �́A�E�F�[�u���b�g�W���ɓK�p����X���b�V���z�[���h�ł��B
%       ���x���ˑ��̎�@�ɑ΂��āATHR_VAL(i,j) �́A���x�� j �ɂ����� 
%       detail �W���ɓK�p����X���b�V���z�[���h�ł��B
%     - L2_Perf(i) �́A���k��ɕۑ������G�l���M�[ (L2-�m����) �̊����ł��B
%     - N0_Perf(i) �́A���k��ɓ����� 0 �̊����ł��B
%
%   ����� 3 �̃I�v�V�������͂��g�p�ł��܂��B
%       [...] = MSWCMPTP(...,S_OR_H) �܂���
%       [...] = MSWCMPTP(...,S_OR_H,KEEPAPP) �܂���
%       [...] = MSWCMPTP(...,S_OR_H,KEEPAPP,IDXSIG)
%       - S_or_H  ('s' �܂��� 'h') �́A�\�t�g�A�܂��̓n�[�h�X���b�V��
%         �z�[���h���Ƃ����Ӗ��ł� (�ڍׂɂ��Ă� MSWTHRESH ���Q��)�B
%       - KEEPAPP (true �܂��� false)�BKEEPAPP �� true �ɓ������ꍇ�A
%         approximation �W���͕ێ�����܂��B
%       - IDXSIG �́A�����M���̃C���f�b�N�X���܂ރx�N�g���A�܂���
%         ������ 'all' �ł��B
%   �f�t�H���g�́A���ꂼ����̂悤�ɂȂ�܂�: 'h', false, 'all'.
%
%   �Q�l mdwtdec, mdwtrec, ddencmp, wdencmp


%   Copyright 1995-2007 The MathWorks, Inc.
