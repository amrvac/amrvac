	.file	"mod_usr.f"
	.text
.Ltext0:
	.p2align 4,,15
	.globl	__mod_usr_MOD_ebr_init
	.type	__mod_usr_MOD_ebr_init, @function
__mod_usr_MOD_ebr_init:
.LVL0:
.LFB0:
	.file 1 "mod_usr.f"
	.loc 1 15 0 view -0
	.cfi_startproc
	.loc 1 15 0 is_stmt 0 view .LVU1
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$184, %rsp
	.cfi_def_cfa_offset 240
	.loc 1 16 0 is_stmt 1 view .LVU2
	movslq	(%rdi), %rdi
.LVL1:
	.loc 1 16 0 is_stmt 0 view .LVU3
	movslq	(%rdx), %r10
.LVL2:
	.loc 1 16 0 view .LVU4
	xorl	%edx, %edx
.LVL3:
	.loc 1 15 0 is_stmt 1 view .LVU5
	movq	256(%rsp), %rax
	.loc 1 16 0 view .LVU6
	movslq	(%rsi), %rsi
.LVL4:
	.loc 1 16 0 is_stmt 0 view .LVU7
	movslq	(%rcx), %rcx
.LVL5:
	.loc 1 23 0 is_stmt 1 view .LVU8
	movslq	(%r8), %rbp
.LVL6:
	.loc 1 15 0 view .LVU9
	movq	%rax, 24(%rsp)
	.loc 1 16 0 view .LVU10
	movl	$1, %eax
	.loc 1 23 0 view .LVU11
	movslq	(%r9), %r15
	.loc 1 16 0 view .LVU12
	movq	%rax, %r14
	.loc 1 23 0 view .LVU13
	movq	%rax, %r12
	.loc 1 16 0 view .LVU14
	subq	%rdi, %r14
	addq	%r10, %r14
	movq	%rax, %r10
.LVL7:
	.loc 1 16 0 is_stmt 0 view .LVU15
	cmovs	%rdx, %r14
	subq	%rsi, %r10
	addq	%rcx, %r10
	.loc 1 23 0 is_stmt 1 view .LVU16
	movq	240(%rsp), %rcx
.LVL8:
	.loc 1 16 0 view .LVU17
	imulq	%r14, %r10
	testq	%r10, %r10
	movq	%r10, %rbx
	cmovs	%rdx, %rbx
	imulq	%r14, %rsi
.LVL9:
	.loc 1 16 0 is_stmt 0 view .LVU18
	negq	%rdi
.LVL10:
	.loc 1 23 0 is_stmt 1 view .LVU19
	subq	%rbp, %r12
	.loc 1 16 0 view .LVU20
	movq	%rbx, 8(%rsp)
	subq	%rsi, %rdi
.LVL11:
	.loc 1 16 0 is_stmt 0 view .LVU21
	subq	%rbx, %rdi
	.loc 1 23 0 is_stmt 1 view .LVU22
	movslq	(%rcx), %rbx
.LVL12:
	.loc 1 23 0 is_stmt 0 view .LVU23
	movq	248(%rsp), %rcx
	.loc 1 16 0 is_stmt 1 view .LVU24
	movq	%rdi, 16(%rsp)
	.loc 1 23 0 view .LVU25
	addq	%rbx, %r12
	movslq	(%rcx), %r13
.LVL13:
	.loc 1 23 0 is_stmt 0 view .LVU26
	movq	%rax, %rcx
	cmovs	%rdx, %r12
	subq	%r15, %rcx
	movq	%rcx, %rdi
	addq	%r13, %rdi
	imulq	%r12, %rdi
	testq	%rdi, %rdi
	cmovs	%rdx, %rdi
	salq	$3, %rdi
	cmove	%rax, %rdi
	call	malloc
.LVL14:
	.loc 1 26 0 is_stmt 1 view .LVU27
	movsd	__mod_ard_phys_MOD_br_a(%rip), %xmm3
	.loc 1 23 0 view .LVU28
	movq	%rbp, %rdx
	.loc 1 27 0 view .LVU29
	movsd	__mod_ard_phys_MOD_br_b(%rip), %xmm2
	.loc 1 23 0 view .LVU30
	movq	%rax, %r8
.LVL15:
	.loc 1 23 0 is_stmt 0 view .LVU31
	movq	%r12, %rax
.LVL16:
	.loc 1 23 0 view .LVU32
	negq	%rdx
.LBB2:
	.loc 1 30 0 is_stmt 1 view .LVU33
	movl	$770, %ecx
.LBE2:
	.loc 1 27 0 view .LVU34
	divsd	%xmm3, %xmm2
	.loc 1 23 0 view .LVU35
	imulq	%r15, %rax
	movq	%rdx, %rsi
.LBB3:
	.loc 1 30 0 view .LVU36
	movq	%r8, 80(%rsp)
.LBE3:
	.loc 1 28 0 view .LVU37
	movapd	%xmm3, %xmm1
.LBB4:
	.loc 1 30 0 view .LVU38
	leaq	80(%rsp), %rdi
	movq	%r8, 40(%rsp)
.LBE4:
	.loc 1 28 0 view .LVU39
	mulsd	__mod_ard_phys_MOD_br_c(%rip), %xmm1
	.loc 1 23 0 view .LVU40
	movq	%rdx, 72(%rsp)
	subq	%rax, %rsi
.LBB5:
	.loc 1 30 0 view .LVU41
	xorl	%eax, %eax
.LBE5:
	.loc 1 28 0 view .LVU42
	divsd	__mod_ard_phys_MOD_br_d(%rip), %xmm1
.LBB6:
	.loc 1 30 0 view .LVU43
	movq	$0, 104(%rsp)
.LBE6:
	.loc 1 28 0 view .LVU44
	movsd	%xmm1, 48(%rsp)
	.loc 1 23 0 view .LVU45
	movq	%rsi, 32(%rsp)
	.loc 1 28 0 view .LVU46
	movsd	%xmm3, 64(%rsp)
.LBB7:
	.loc 1 30 0 view .LVU47
	movq	$8, 112(%rsp)
	movq	$8, 96(%rsp)
	movw	%cx, 108(%rsp)
	movq	%rbp, 128(%rsp)
	movq	%rbx, 136(%rsp)
	movq	$1, 120(%rsp)
	movq	%r15, 152(%rsp)
	movq	%r13, 160(%rsp)
	movq	%r12, 144(%rsp)
	movq	%rsi, 88(%rsp)
.LBE7:
	.loc 1 27 0 view .LVU48
	movsd	%xmm2, 56(%rsp)
.LVL17:
.LBB8:
	.loc 1 30 0 view .LVU49
	call	_gfortran_arandom_r8
.LVL18:
	.loc 1 30 0 is_stmt 0 view .LVU50
.LBE8:
.LBB9:
.LBB10:
	.loc 1 31 0 is_stmt 1 view .LVU51
	movq	40(%rsp), %r8
.LBE10:
	movslq	__mod_ard_phys_MOD_u_(%rip), %rax
.LBB12:
	imulq	8(%rsp), %rax
	addq	16(%rsp), %rax
	cmpq	%r13, %r15
	movsd	48(%rsp), %xmm1
	movsd	56(%rsp), %xmm2
	jg	.L2
	movq	72(%rsp), %rdx
	movq	24(%rsp), %rcx
	leaq	1(%r15), %rdi
	leaq	0(,%r12,8), %r11
	movsd	.LC0(%rip), %xmm4
	movsd	64(%rsp), %xmm3
	leaq	0(,%r14,8), %r10
	leaq	2(%r13), %r9
	leaq	(%r8,%rdx,8), %rsi
	movq	%r14, %rdx
	imulq	%r15, %rdx
	addq	%rdx, %rax
	leaq	(%rcx,%rax,8), %rcx
	.p2align 4,,10
	.p2align 3
.L5:
.LBB11:
	.loc 1 31 0 is_stmt 0 discriminator 2 view .LVU52
	cmpq	%rbx, %rbp
	jg	.L3
	movq	%rbp, %rax
	jmp	.L4
	.p2align 4,,10
	.p2align 3
.L14:
	.loc 1 31 0 discriminator 2 view .LVU53
	movq	%rdx, %rax
.L4:
	.loc 1 31 0 discriminator 5 view .LVU54
	movsd	(%rsi,%rax,8), %xmm0
	leaq	1(%rax), %rdx
	mulsd	%xmm4, %xmm0
	addsd	%xmm3, %xmm0
	movsd	%xmm0, (%rcx,%rax,8)
	cmpq	%rbx, %rax
	jne	.L14
.L3:
	.loc 1 31 0 discriminator 5 view .LVU55
	addq	$1, %rdi
	addq	%r11, %rsi
	addq	%r10, %rcx
.LBE11:
	.loc 1 31 0 view .LVU56
	cmpq	%r9, %rdi
	jne	.L5
.L2:
	.loc 1 31 0 view .LVU57
.LBE12:
.LBE9:
.LBB13:
	.loc 1 32 0 is_stmt 1 view .LVU58
	movq	32(%rsp), %rax
	movl	$770, %edx
	leaq	80(%rsp), %rdi
	movsd	%xmm1, 48(%rsp)
.LVL19:
	.loc 1 32 0 is_stmt 0 view .LVU59
	movq	$0, 104(%rsp)
	movq	%rax, 88(%rsp)
	xorl	%eax, %eax
	movq	%r8, 80(%rsp)
	movq	%r8, 40(%rsp)
.LVL20:
	.loc 1 32 0 view .LVU60
	movsd	%xmm2, 56(%rsp)
.LVL21:
	.loc 1 32 0 view .LVU61
	movq	$8, 112(%rsp)
	movq	$8, 96(%rsp)
	movw	%dx, 108(%rsp)
	movq	%rbp, 128(%rsp)
	movq	%rbx, 136(%rsp)
	movq	$1, 120(%rsp)
	movq	%r15, 152(%rsp)
	movq	%r13, 160(%rsp)
	movq	%r12, 144(%rsp)
	call	_gfortran_arandom_r8
.LVL22:
.LBE13:
.LBB14:
.LBB15:
	.loc 1 33 0 is_stmt 1 view .LVU62
	movq	40(%rsp), %r8
.LBE15:
	movslq	__mod_ard_phys_MOD_v_(%rip), %rax
.LBB17:
	imulq	8(%rsp), %rax
	addq	16(%rsp), %rax
	cmpq	%r13, %r15
	movsd	48(%rsp), %xmm1
	jg	.L6
	movq	%rbp, %rdx
	movq	24(%rsp), %rcx
	leaq	1(%r15), %rdi
	leaq	2(%r13), %r9
	negq	%rdx
	movsd	.LC0(%rip), %xmm3
	movsd	56(%rsp), %xmm2
	leaq	0(,%r12,8), %r11
	leaq	(%r8,%rdx,8), %rsi
	movq	%r14, %rdx
	leaq	0(,%r14,8), %r10
	imulq	%r15, %rdx
	addq	%rdx, %rax
	leaq	(%rcx,%rax,8), %rcx
	.p2align 4,,10
	.p2align 3
.L9:
.LBB16:
	.loc 1 33 0 is_stmt 0 discriminator 2 view .LVU63
	cmpq	%rbx, %rbp
	jg	.L7
	movq	%rbp, %rax
	jmp	.L8
	.p2align 4,,10
	.p2align 3
.L15:
	.loc 1 33 0 discriminator 2 view .LVU64
	movq	%rdx, %rax
.L8:
	.loc 1 33 0 discriminator 5 view .LVU65
	movsd	(%rsi,%rax,8), %xmm0
	leaq	1(%rax), %rdx
	mulsd	%xmm3, %xmm0
	addsd	%xmm2, %xmm0
	movsd	%xmm0, (%rcx,%rax,8)
	cmpq	%rbx, %rax
	jne	.L15
.L7:
	.loc 1 33 0 discriminator 5 view .LVU66
	addq	$1, %rdi
	addq	%r11, %rsi
	addq	%r10, %rcx
.LBE16:
	.loc 1 33 0 view .LVU67
	cmpq	%r9, %rdi
	jne	.L9
.L6:
	.loc 1 33 0 view .LVU68
.LBE17:
.LBE14:
.LBB18:
	.loc 1 34 0 is_stmt 1 view .LVU69
	movl	$770, %eax
	leaq	80(%rsp), %rdi
	movq	%r8, 80(%rsp)
	movq	$0, 104(%rsp)
	movw	%ax, 108(%rsp)
	movq	32(%rsp), %rax
	movq	%r8, 40(%rsp)
	movq	%rax, 88(%rsp)
	xorl	%eax, %eax
	movsd	%xmm1, 48(%rsp)
	movq	$8, 112(%rsp)
	movq	$8, 96(%rsp)
	movq	%rbp, 128(%rsp)
	movq	%rbx, 136(%rsp)
	movq	$1, 120(%rsp)
	movq	%r15, 152(%rsp)
	movq	%r13, 160(%rsp)
	movq	%r12, 144(%rsp)
	call	_gfortran_arandom_r8
.LVL23:
.LBE18:
.LBB19:
	.loc 1 35 0 view .LVU70
	movslq	__mod_ard_phys_MOD_w_(%rip), %rax
.LBB20:
	imulq	8(%rsp), %rax
	addq	16(%rsp), %rax
	cmpq	%r13, %r15
	movq	40(%rsp), %r8
	jg	.L10
	movq	%rbp, %rdx
	movq	24(%rsp), %rcx
	leaq	1(%r15), %rdi
	salq	$3, %r12
	leaq	0(,%r14,8), %r9
	negq	%rdx
	movsd	48(%rsp), %xmm1
	movsd	.LC0(%rip), %xmm2
	imulq	%r15, %r14
	leaq	(%r8,%rdx,8), %rsi
	addq	$2, %r13
.LVL24:
	.loc 1 35 0 is_stmt 0 view .LVU71
	addq	%r14, %rax
	leaq	(%rcx,%rax,8), %rcx
	.p2align 4,,10
	.p2align 3
.L13:
.LBB21:
	.loc 1 35 0 discriminator 2 view .LVU72
	cmpq	%rbx, %rbp
	jg	.L11
	movq	%rbp, %rax
	jmp	.L12
	.p2align 4,,10
	.p2align 3
.L16:
	.loc 1 35 0 discriminator 2 view .LVU73
	movq	%rdx, %rax
.L12:
	.loc 1 35 0 discriminator 5 view .LVU74
	movsd	(%rsi,%rax,8), %xmm0
	leaq	1(%rax), %rdx
	mulsd	%xmm2, %xmm0
	addsd	%xmm1, %xmm0
	movsd	%xmm0, (%rcx,%rax,8)
	cmpq	%rax, %rbx
	jne	.L16
.L11:
	.loc 1 35 0 discriminator 5 view .LVU75
	addq	$1, %rdi
	addq	%r12, %rsi
	addq	%r9, %rcx
.LBE21:
	.loc 1 35 0 view .LVU76
	cmpq	%r13, %rdi
	jne	.L13
.LVL25:
.L10:
	.loc 1 35 0 view .LVU77
.LBE20:
.LBE19:
	.loc 1 36 0 is_stmt 1 view .LVU78
	addq	$184, %rsp
	.cfi_def_cfa_offset 56
	.loc 1 23 0 view .LVU79
	movq	%r8, %rdi
	.loc 1 36 0 view .LVU80
	popq	%rbx
	.cfi_def_cfa_offset 48
.LVL26:
	.loc 1 36 0 is_stmt 0 view .LVU81
	popq	%rbp
	.cfi_def_cfa_offset 40
.LVL27:
	.loc 1 36 0 view .LVU82
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
.LVL28:
	.loc 1 23 0 is_stmt 1 view .LVU83
	jmp	free
.LVL29:
	.cfi_endproc
.LFE0:
	.size	__mod_usr_MOD_ebr_init, .-__mod_usr_MOD_ebr_init
	.section	.rodata
.LC1:
	.ascii	"Cartesian"
	.text
	.p2align 4,,15
	.globl	__mod_usr_MOD_usr_init
	.type	__mod_usr_MOD_usr_init, @function
__mod_usr_MOD_usr_init:
.LFB1:
	.loc 1 8 0 view -0
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	.loc 1 9 0 view .LVU85
	movl	$9, %esi
	movl	$.LC1, %edi
	call	__mod_geometry_MOD_set_coordinate_system
.LVL30:
	.loc 1 10 0 view .LVU86
	call	__mod_ard_MOD_ard_activate
.LVL31:
	.loc 1 12 0 view .LVU87
	movq	$__mod_usr_MOD_ebr_init, __mod_usr_methods_MOD_usr_init_one_grid(%rip)
	.loc 1 13 0 view .LVU88
	addq	$8, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE1:
	.size	__mod_usr_MOD_usr_init, .-__mod_usr_MOD_usr_init
	.comm	mpi_fortran_weights_empty_,4,16
	.comm	mpi_fortran_unweighted_,4,16
	.comm	mpi_fortran_statuses_ignore_,24,16
	.comm	mpi_fortran_status_ignore_,24,16
	.comm	mpi_fortran_in_place_,4,16
	.comm	mpi_fortran_errcodes_ignore_,4,16
	.comm	mpi_fortran_bottom_,4,16
	.comm	mpi_fortran_argvs_null_,1,16
	.comm	mpi_fortran_argv_null_,1,16
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC0:
	.long	2576980378
	.long	1069128089
	.text
.Letext0:
	.file 2 "<built-in>"
	.section	.debug_info,"",@progbits
.Ldebug_info0:
	.long	0x3a2
	.value	0x4
	.long	.Ldebug_abbrev0
	.byte	0x8
	.uleb128 0x1
	.long	.LASF21
	.byte	0xe
	.byte	0x2
	.long	.LASF22
	.long	.LASF23
	.quad	.Ltext0
	.quad	.Letext0-.Ltext0
	.long	.Ldebug_line0
	.uleb128 0x2
	.long	.LASF24
	.byte	0x1
	.byte	0x1
	.long	0x2f1
	.uleb128 0x3
	.byte	0x1
	.byte	0x3
	.long	0x2f1
	.uleb128 0x4
	.long	.LASF25
	.byte	0x1
	.byte	0x8
	.long	.LASF26
	.quad	.LFB1
	.quad	.LFE1-.LFB1
	.uleb128 0x1
	.byte	0x9c
	.long	0x93
	.uleb128 0x5
	.quad	.LVL30
	.long	0x36e
	.long	0x85
	.uleb128 0x6
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x9
	.byte	0x3
	.quad	.LC1
	.uleb128 0x6
	.uleb128 0x1
	.byte	0x54
	.uleb128 0x1
	.byte	0x39
	.byte	0
	.uleb128 0x7
	.quad	.LVL31
	.long	0x379
	.byte	0
	.uleb128 0x8
	.long	.LASF27
	.byte	0x1
	.byte	0xf
	.long	.LASF28
	.quad	.LFB0
	.quad	.LFE0-.LFB0
	.uleb128 0x1
	.byte	0x9c
	.uleb128 0x9
	.long	.LASF0
	.byte	0x1
	.byte	0xf
	.long	0x2fd
	.long	.LLST0
	.long	.LVUS0
	.uleb128 0x9
	.long	.LASF1
	.byte	0x1
	.byte	0xf
	.long	0x2fd
	.long	.LLST1
	.long	.LVUS1
	.uleb128 0x9
	.long	.LASF2
	.byte	0x1
	.byte	0xf
	.long	0x2fd
	.long	.LLST2
	.long	.LVUS2
	.uleb128 0x9
	.long	.LASF3
	.byte	0x1
	.byte	0xf
	.long	0x2fd
	.long	.LLST3
	.long	.LVUS3
	.uleb128 0x9
	.long	.LASF4
	.byte	0x1
	.byte	0xf
	.long	0x2fd
	.long	.LLST4
	.long	.LVUS4
	.uleb128 0x9
	.long	.LASF5
	.byte	0x1
	.byte	0xf
	.long	0x2fd
	.long	.LLST5
	.long	.LVUS5
	.uleb128 0xa
	.long	.LASF6
	.byte	0x1
	.byte	0xf
	.long	0x2fd
	.uleb128 0x3
	.byte	0x91
	.sleb128 0
	.byte	0x6
	.uleb128 0xa
	.long	.LASF7
	.byte	0x1
	.byte	0xf
	.long	0x2fd
	.uleb128 0x3
	.byte	0x91
	.sleb128 8
	.byte	0x6
	.uleb128 0xb
	.string	"w"
	.byte	0x1
	.byte	0xf
	.long	0x304
	.uleb128 0x3
	.byte	0x91
	.sleb128 16
	.byte	0x6
	.uleb128 0xb
	.string	"x"
	.byte	0x1
	.byte	0xf
	.long	0x339
	.uleb128 0x3
	.byte	0x91
	.sleb128 24
	.byte	0x6
	.uleb128 0xc
	.string	"ssu"
	.byte	0x1
	.byte	0x18
	.long	0x332
	.long	.LLST6
	.long	.LVUS6
	.uleb128 0xc
	.string	"ssv"
	.byte	0x1
	.byte	0x18
	.long	0x332
	.long	.LLST7
	.long	.LVUS7
	.uleb128 0xc
	.string	"ssw"
	.byte	0x1
	.byte	0x18
	.long	0x332
	.long	.LLST8
	.long	.LVUS8
	.uleb128 0xd
	.long	0x2f6
	.long	.LLST9
	.long	.LVUS9
	.uleb128 0xd
	.long	0x2f6
	.long	.LLST10
	.long	.LVUS10
	.uleb128 0xd
	.long	0x2f6
	.long	.LLST11
	.long	.LVUS11
	.uleb128 0xd
	.long	0x2f6
	.long	.LLST12
	.long	.LVUS12
	.uleb128 0xe
	.long	.LASF8
	.byte	0x1
	.byte	0x17
	.long	0x349
	.long	.LLST13
	.long	.LVUS13
	.uleb128 0xd
	.long	0x2f6
	.long	.LLST14
	.long	.LVUS14
	.uleb128 0xd
	.long	0x2f6
	.long	.LLST15
	.long	.LVUS15
	.uleb128 0xd
	.long	0x2f6
	.long	.LLST16
	.long	.LVUS16
	.uleb128 0xd
	.long	0x2f6
	.long	.LLST17
	.long	.LVUS17
	.uleb128 0xf
	.long	0x2f6
	.uleb128 0x10
	.long	0x22e
	.uleb128 0x11
	.quad	.LVL18
	.long	0x384
	.uleb128 0x6
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x3
	.byte	0x91
	.sleb128 -160
	.byte	0
	.byte	0
	.uleb128 0x10
	.long	0x237
	.uleb128 0x12
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x10
	.long	0x252
	.uleb128 0x11
	.quad	.LVL22
	.long	0x384
	.uleb128 0x6
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x3
	.byte	0x91
	.sleb128 -160
	.byte	0
	.byte	0
	.uleb128 0x10
	.long	0x25b
	.uleb128 0x12
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x10
	.long	0x276
	.uleb128 0x11
	.quad	.LVL23
	.long	0x384
	.uleb128 0x6
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x3
	.byte	0x91
	.sleb128 -160
	.byte	0
	.byte	0
	.uleb128 0x10
	.long	0x27f
	.uleb128 0x12
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x5
	.quad	.LVL14
	.long	0x38f
	.long	0x2d9
	.uleb128 0x6
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x44
	.byte	0x7d
	.sleb128 0
	.byte	0x7f
	.sleb128 0
	.byte	0x1c
	.byte	0x23
	.uleb128 0x1
	.byte	0x7c
	.sleb128 0
	.byte	0x1e
	.byte	0x30
	.byte	0x7d
	.sleb128 0
	.byte	0x7f
	.sleb128 0
	.byte	0x1c
	.byte	0x23
	.uleb128 0x1
	.byte	0x7c
	.sleb128 0
	.byte	0x1e
	.byte	0x30
	.byte	0x2a
	.byte	0x28
	.value	0x1
	.byte	0x16
	.byte	0x13
	.byte	0x33
	.byte	0x24
	.byte	0x31
	.byte	0x7d
	.sleb128 0
	.byte	0x7f
	.sleb128 0
	.byte	0x1c
	.byte	0x23
	.uleb128 0x1
	.byte	0x7c
	.sleb128 0
	.byte	0x1e
	.byte	0x30
	.byte	0x7d
	.sleb128 0
	.byte	0x7f
	.sleb128 0
	.byte	0x1c
	.byte	0x23
	.uleb128 0x1
	.byte	0x7c
	.sleb128 0
	.byte	0x1e
	.byte	0x30
	.byte	0x2a
	.byte	0x28
	.value	0x1
	.byte	0x16
	.byte	0x13
	.byte	0x33
	.byte	0x24
	.byte	0x30
	.byte	0x2e
	.byte	0x28
	.value	0x1
	.byte	0x16
	.byte	0x13
	.byte	0
	.uleb128 0x14
	.quad	.LVL29
	.long	0x39a
	.uleb128 0x6
	.uleb128 0x1
	.byte	0x55
	.uleb128 0x4
	.byte	0x91
	.sleb128 -200
	.byte	0x6
	.byte	0
	.byte	0
	.byte	0
	.uleb128 0x15
	.long	.LASF29
	.uleb128 0x16
	.byte	0x8
	.byte	0x5
	.long	.LASF9
	.uleb128 0x16
	.byte	0x4
	.byte	0x5
	.long	.LASF10
	.uleb128 0x17
	.byte	0x1
	.long	0x332
	.long	0x332
	.uleb128 0x18
	.long	0x2f6
	.long	0x1da
	.long	0x1e7
	.uleb128 0x18
	.long	0x2f6
	.long	0x1f4
	.long	0x201
	.uleb128 0x19
	.long	0x2f6
	.long	0x20e
	.byte	0
	.uleb128 0x16
	.byte	0x8
	.byte	0x4
	.long	.LASF11
	.uleb128 0x1a
	.long	0x332
	.long	0x349
	.uleb128 0x1b
	.long	0x2f6
	.sleb128 0
	.byte	0
	.uleb128 0x17
	.byte	0x1
	.long	0x332
	.long	0x36e
	.uleb128 0x18
	.long	0x2f6
	.long	0x193
	.long	0x1a0
	.uleb128 0x18
	.long	0x2f6
	.long	0x1ad
	.long	0x1ba
	.byte	0
	.uleb128 0x1c
	.long	.LASF12
	.long	.LASF14
	.byte	0x1
	.byte	0x9
	.uleb128 0x1c
	.long	.LASF13
	.long	.LASF15
	.byte	0x1
	.byte	0xa
	.uleb128 0x1c
	.long	.LASF16
	.long	.LASF16
	.byte	0x1
	.byte	0x1e
	.uleb128 0x1c
	.long	.LASF17
	.long	.LASF18
	.byte	0x2
	.byte	0
	.uleb128 0x1c
	.long	.LASF19
	.long	.LASF20
	.byte	0x2
	.byte	0
	.byte	0
	.section	.debug_abbrev,"",@progbits
.Ldebug_abbrev0:
	.uleb128 0x1
	.uleb128 0x11
	.byte	0x1
	.uleb128 0x25
	.uleb128 0xe
	.uleb128 0x13
	.uleb128 0xb
	.uleb128 0x42
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x1b
	.uleb128 0xe
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x10
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0x2
	.uleb128 0x1e
	.byte	0x1
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x3
	.uleb128 0x3a
	.byte	0
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x18
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x4
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x6e
	.uleb128 0xe
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x40
	.uleb128 0x18
	.uleb128 0x2117
	.uleb128 0x19
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x5
	.uleb128 0x4109
	.byte	0x1
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x31
	.uleb128 0x13
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x6
	.uleb128 0x410a
	.byte	0
	.uleb128 0x2
	.uleb128 0x18
	.uleb128 0x2111
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0x7
	.uleb128 0x4109
	.byte	0
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x31
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x8
	.uleb128 0x2e
	.byte	0x1
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x6e
	.uleb128 0xe
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x12
	.uleb128 0x7
	.uleb128 0x40
	.uleb128 0x18
	.uleb128 0x2117
	.uleb128 0x19
	.byte	0
	.byte	0
	.uleb128 0x9
	.uleb128 0x5
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0xa
	.uleb128 0x5
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0xb
	.uleb128 0x5
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x18
	.byte	0
	.byte	0
	.uleb128 0xc
	.uleb128 0x34
	.byte	0
	.uleb128 0x3
	.uleb128 0x8
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0xd
	.uleb128 0x34
	.byte	0
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x34
	.uleb128 0x19
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0xe
	.uleb128 0x34
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2
	.uleb128 0x17
	.uleb128 0x2137
	.uleb128 0x17
	.byte	0
	.byte	0
	.uleb128 0xf
	.uleb128 0x34
	.byte	0
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x34
	.uleb128 0x19
	.byte	0
	.byte	0
	.uleb128 0x10
	.uleb128 0xb
	.byte	0x1
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x11
	.uleb128 0x4109
	.byte	0x1
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x31
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x12
	.uleb128 0xb
	.byte	0x1
	.byte	0
	.byte	0
	.uleb128 0x13
	.uleb128 0xb
	.byte	0
	.byte	0
	.byte	0
	.uleb128 0x14
	.uleb128 0x4109
	.byte	0x1
	.uleb128 0x11
	.uleb128 0x1
	.uleb128 0x2115
	.uleb128 0x19
	.uleb128 0x31
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x15
	.uleb128 0x1e
	.byte	0
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3c
	.uleb128 0x19
	.byte	0
	.byte	0
	.uleb128 0x16
	.uleb128 0x24
	.byte	0
	.uleb128 0xb
	.uleb128 0xb
	.uleb128 0x3e
	.uleb128 0xb
	.uleb128 0x3
	.uleb128 0xe
	.byte	0
	.byte	0
	.uleb128 0x17
	.uleb128 0x1
	.byte	0x1
	.uleb128 0x9
	.uleb128 0xb
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x18
	.uleb128 0x21
	.byte	0
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x22
	.uleb128 0x13
	.uleb128 0x2f
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x19
	.uleb128 0x21
	.byte	0
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x2f
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x1a
	.uleb128 0x1
	.byte	0x1
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x1
	.uleb128 0x13
	.byte	0
	.byte	0
	.uleb128 0x1b
	.uleb128 0x21
	.byte	0
	.uleb128 0x49
	.uleb128 0x13
	.uleb128 0x22
	.uleb128 0xd
	.byte	0
	.byte	0
	.uleb128 0x1c
	.uleb128 0x2e
	.byte	0
	.uleb128 0x3f
	.uleb128 0x19
	.uleb128 0x3c
	.uleb128 0x19
	.uleb128 0x6e
	.uleb128 0xe
	.uleb128 0x3
	.uleb128 0xe
	.uleb128 0x3a
	.uleb128 0xb
	.uleb128 0x3b
	.uleb128 0xb
	.byte	0
	.byte	0
	.byte	0
	.section	.debug_loc,"",@progbits
.Ldebug_loc0:
.LVUS0:
	.uleb128 0
	.uleb128 .LVU3
	.uleb128 .LVU3
	.uleb128 0
.LLST0:
	.quad	.LVL0-.Ltext0
	.quad	.LVL1-.Ltext0
	.value	0x2
	.byte	0x75
	.sleb128 0
	.quad	.LVL1-.Ltext0
	.quad	.LFE0-.Ltext0
	.value	0x3
	.byte	0xf3
	.uleb128 0x1
	.byte	0x55
	.quad	0
	.quad	0
.LVUS1:
	.uleb128 0
	.uleb128 .LVU7
	.uleb128 .LVU7
	.uleb128 0
.LLST1:
	.quad	.LVL0-.Ltext0
	.quad	.LVL4-.Ltext0
	.value	0x2
	.byte	0x74
	.sleb128 0
	.quad	.LVL4-.Ltext0
	.quad	.LFE0-.Ltext0
	.value	0x3
	.byte	0xf3
	.uleb128 0x1
	.byte	0x54
	.quad	0
	.quad	0
.LVUS2:
	.uleb128 0
	.uleb128 .LVU5
	.uleb128 .LVU5
	.uleb128 0
.LLST2:
	.quad	.LVL0-.Ltext0
	.quad	.LVL3-.Ltext0
	.value	0x2
	.byte	0x71
	.sleb128 0
	.quad	.LVL3-.Ltext0
	.quad	.LFE0-.Ltext0
	.value	0x3
	.byte	0xf3
	.uleb128 0x1
	.byte	0x51
	.quad	0
	.quad	0
.LVUS3:
	.uleb128 0
	.uleb128 .LVU8
	.uleb128 .LVU8
	.uleb128 0
.LLST3:
	.quad	.LVL0-.Ltext0
	.quad	.LVL5-.Ltext0
	.value	0x2
	.byte	0x72
	.sleb128 0
	.quad	.LVL5-.Ltext0
	.quad	.LFE0-.Ltext0
	.value	0x3
	.byte	0xf3
	.uleb128 0x1
	.byte	0x52
	.quad	0
	.quad	0
.LVUS4:
	.uleb128 0
	.uleb128 .LVU27
	.uleb128 .LVU27
	.uleb128 0
.LLST4:
	.quad	.LVL0-.Ltext0
	.quad	.LVL14-1-.Ltext0
	.value	0x2
	.byte	0x78
	.sleb128 0
	.quad	.LVL14-1-.Ltext0
	.quad	.LFE0-.Ltext0
	.value	0x3
	.byte	0xf3
	.uleb128 0x1
	.byte	0x58
	.quad	0
	.quad	0
.LVUS5:
	.uleb128 0
	.uleb128 .LVU27
	.uleb128 .LVU27
	.uleb128 0
.LLST5:
	.quad	.LVL0-.Ltext0
	.quad	.LVL14-1-.Ltext0
	.value	0x2
	.byte	0x79
	.sleb128 0
	.quad	.LVL14-1-.Ltext0
	.quad	.LFE0-.Ltext0
	.value	0x3
	.byte	0xf3
	.uleb128 0x1
	.byte	0x59
	.quad	0
	.quad	0
.LVUS6:
	.uleb128 .LVU31
	.uleb128 .LVU50
	.uleb128 .LVU50
	.uleb128 0
.LLST6:
	.quad	.LVL15-.Ltext0
	.quad	.LVL18-1-.Ltext0
	.value	0x1
	.byte	0x64
	.quad	.LVL18-1-.Ltext0
	.quad	.LFE0-.Ltext0
	.value	0x3
	.byte	0x91
	.sleb128 -176
	.quad	0
	.quad	0
.LVUS7:
	.uleb128 .LVU49
	.uleb128 .LVU50
	.uleb128 .LVU50
	.uleb128 .LVU61
.LLST7:
	.quad	.LVL17-.Ltext0
	.quad	.LVL18-1-.Ltext0
	.value	0x1
	.byte	0x63
	.quad	.LVL18-1-.Ltext0
	.quad	.LVL21-.Ltext0
	.value	0x3
	.byte	0x91
	.sleb128 -184
	.quad	0
	.quad	0
.LVUS8:
	.uleb128 .LVU49
	.uleb128 .LVU50
	.uleb128 .LVU50
	.uleb128 .LVU59
.LLST8:
	.quad	.LVL17-.Ltext0
	.quad	.LVL18-1-.Ltext0
	.value	0x1
	.byte	0x62
	.quad	.LVL18-1-.Ltext0
	.quad	.LVL19-.Ltext0
	.value	0x3
	.byte	0x91
	.sleb128 -192
	.quad	0
	.quad	0
.LVUS9:
	.uleb128 .LVU9
	.uleb128 .LVU82
	.uleb128 .LVU82
	.uleb128 0
.LLST9:
	.quad	.LVL6-.Ltext0
	.quad	.LVL27-.Ltext0
	.value	0x1
	.byte	0x56
	.quad	.LVL27-.Ltext0
	.quad	.LFE0-.Ltext0
	.value	0x6
	.byte	0x91
	.sleb128 -168
	.byte	0x6
	.byte	0x1f
	.byte	0x9f
	.quad	0
	.quad	0
.LVUS10:
	.uleb128 .LVU23
	.uleb128 .LVU81
.LLST10:
	.quad	.LVL12-.Ltext0
	.quad	.LVL26-.Ltext0
	.value	0x1
	.byte	0x53
	.quad	0
	.quad	0
.LVUS11:
	.uleb128 .LVU23
	.uleb128 .LVU83
.LLST11:
	.quad	.LVL12-.Ltext0
	.quad	.LVL28-.Ltext0
	.value	0x1
	.byte	0x5f
	.quad	0
	.quad	0
.LVUS12:
	.uleb128 .LVU26
	.uleb128 .LVU71
	.uleb128 .LVU71
	.uleb128 .LVU77
.LLST12:
	.quad	.LVL13-.Ltext0
	.quad	.LVL24-.Ltext0
	.value	0x1
	.byte	0x5d
	.quad	.LVL24-.Ltext0
	.quad	.LVL25-.Ltext0
	.value	0x3
	.byte	0x7d
	.sleb128 -2
	.byte	0x9f
	.quad	0
	.quad	0
.LVUS13:
	.uleb128 .LVU31
	.uleb128 .LVU32
	.uleb128 .LVU32
	.uleb128 .LVU50
	.uleb128 .LVU50
	.uleb128 .LVU60
.LLST13:
	.quad	.LVL15-.Ltext0
	.quad	.LVL16-.Ltext0
	.value	0x2
	.byte	0x70
	.sleb128 0
	.quad	.LVL16-.Ltext0
	.quad	.LVL18-1-.Ltext0
	.value	0x2
	.byte	0x78
	.sleb128 0
	.quad	.LVL18-1-.Ltext0
	.quad	.LVL20-.Ltext0
	.value	0x4
	.byte	0x91
	.sleb128 -200
	.byte	0x6
	.quad	0
	.quad	0
.LVUS14:
	.uleb128 .LVU3
	.uleb128 .LVU19
	.uleb128 .LVU19
	.uleb128 .LVU21
	.uleb128 .LVU21
	.uleb128 .LVU27
.LLST14:
	.quad	.LVL1-.Ltext0
	.quad	.LVL10-.Ltext0
	.value	0x1
	.byte	0x55
	.quad	.LVL10-.Ltext0
	.quad	.LVL11-.Ltext0
	.value	0x4
	.byte	0x75
	.sleb128 0
	.byte	0x1f
	.byte	0x9f
	.quad	.LVL11-.Ltext0
	.quad	.LVL14-1-.Ltext0
	.value	0xc
	.byte	0xf3
	.uleb128 0x1
	.byte	0x55
	.byte	0x94
	.byte	0x4
	.byte	0x8
	.byte	0x20
	.byte	0x24
	.byte	0x8
	.byte	0x20
	.byte	0x26
	.byte	0x9f
	.quad	0
	.quad	0
.LVUS15:
	.uleb128 .LVU3
	.uleb128 .LVU4
	.uleb128 .LVU4
	.uleb128 .LVU15
	.uleb128 .LVU15
	.uleb128 .LVU27
.LLST15:
	.quad	.LVL1-.Ltext0
	.quad	.LVL2-.Ltext0
	.value	0xb
	.byte	0x71
	.sleb128 0
	.byte	0x94
	.byte	0x4
	.byte	0x8
	.byte	0x20
	.byte	0x24
	.byte	0x8
	.byte	0x20
	.byte	0x26
	.byte	0x9f
	.quad	.LVL2-.Ltext0
	.quad	.LVL7-.Ltext0
	.value	0x1
	.byte	0x5a
	.quad	.LVL7-.Ltext0
	.quad	.LVL14-1-.Ltext0
	.value	0xc
	.byte	0xf3
	.uleb128 0x1
	.byte	0x51
	.byte	0x94
	.byte	0x4
	.byte	0x8
	.byte	0x20
	.byte	0x24
	.byte	0x8
	.byte	0x20
	.byte	0x26
	.byte	0x9f
	.quad	0
	.quad	0
.LVUS16:
	.uleb128 .LVU7
	.uleb128 .LVU18
	.uleb128 .LVU18
	.uleb128 .LVU27
.LLST16:
	.quad	.LVL4-.Ltext0
	.quad	.LVL9-.Ltext0
	.value	0x1
	.byte	0x54
	.quad	.LVL9-.Ltext0
	.quad	.LVL14-1-.Ltext0
	.value	0xc
	.byte	0xf3
	.uleb128 0x1
	.byte	0x54
	.byte	0x94
	.byte	0x4
	.byte	0x8
	.byte	0x20
	.byte	0x24
	.byte	0x8
	.byte	0x20
	.byte	0x26
	.byte	0x9f
	.quad	0
	.quad	0
.LVUS17:
	.uleb128 .LVU7
	.uleb128 .LVU8
	.uleb128 .LVU8
	.uleb128 .LVU17
	.uleb128 .LVU17
	.uleb128 .LVU27
.LLST17:
	.quad	.LVL4-.Ltext0
	.quad	.LVL5-.Ltext0
	.value	0xb
	.byte	0x72
	.sleb128 0
	.byte	0x94
	.byte	0x4
	.byte	0x8
	.byte	0x20
	.byte	0x24
	.byte	0x8
	.byte	0x20
	.byte	0x26
	.byte	0x9f
	.quad	.LVL5-.Ltext0
	.quad	.LVL8-.Ltext0
	.value	0x1
	.byte	0x52
	.quad	.LVL8-.Ltext0
	.quad	.LVL14-1-.Ltext0
	.value	0xc
	.byte	0xf3
	.uleb128 0x1
	.byte	0x52
	.byte	0x94
	.byte	0x4
	.byte	0x8
	.byte	0x20
	.byte	0x24
	.byte	0x8
	.byte	0x20
	.byte	0x26
	.byte	0x9f
	.quad	0
	.quad	0
	.section	.debug_aranges,"",@progbits
	.long	0x2c
	.value	0x2
	.long	.Ldebug_info0
	.byte	0x8
	.byte	0
	.value	0
	.value	0
	.quad	.Ltext0
	.quad	.Letext0-.Ltext0
	.quad	0
	.quad	0
	.section	.debug_line,"",@progbits
.Ldebug_line0:
	.section	.debug_str,"MS",@progbits,1
.LASF21:
	.string	"GNU Fortran2008 8.5.0 20210514 (Red Hat 8.5.0-10) -mtune=generic -march=x86-64 -g -O2 -ffree-form -fintrinsic-modules-path /usr/lib/gcc/x86_64-redhat-linux/8/finclude"
.LASF28:
	.string	"__mod_usr_MOD_ebr_init"
.LASF16:
	.string	"_gfortran_arandom_r8"
.LASF2:
	.string	"ixgmax1"
.LASF3:
	.string	"ixgmax2"
.LASF22:
	.string	"mod_usr.f"
.LASF17:
	.string	"malloc"
.LASF25:
	.string	"usr_init"
.LASF18:
	.string	"__builtin_malloc"
.LASF23:
	.string	"/users/cpa/beatricp/amrvac/tests/ard/ext_brusselator_2d"
.LASF11:
	.string	"real(kind=8)"
.LASF20:
	.string	"__builtin_free"
.LASF26:
	.string	"__mod_usr_MOD_usr_init"
.LASF10:
	.string	"integer(kind=4)"
.LASF13:
	.string	"__mod_ard_MOD_ard_activate"
.LASF4:
	.string	"ixmin1"
.LASF5:
	.string	"ixmin2"
.LASF12:
	.string	"__mod_geometry_MOD_set_coordinate_system"
.LASF8:
	.string	"urand"
.LASF14:
	.string	"set_coordinate_system"
.LASF29:
	.string	"mod_ard"
.LASF24:
	.string	"mod_usr"
.LASF6:
	.string	"ixmax1"
.LASF7:
	.string	"ixmax2"
.LASF9:
	.string	"integer(kind=8)"
.LASF19:
	.string	"free"
.LASF0:
	.string	"ixgmin1"
.LASF1:
	.string	"ixgmin2"
.LASF27:
	.string	"ebr_init"
.LASF15:
	.string	"ard_activate"
	.ident	"GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-10)"
	.section	.note.GNU-stack,"",@progbits
