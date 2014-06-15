/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     PTR_OP = 258,
     INC_OP = 259,
     DEC_OP = 260,
     LEFT_OP = 261,
     RIGHT_OP = 262,
     LE_OP = 263,
     GE_OP = 264,
     EQ_OP = 265,
     NE_OP = 266,
     AND_OP = 267,
     OR_OP = 268,
     MUL_ASSIGN = 269,
     DIV_ASSIGN = 270,
     MOD_ASSIGN = 271,
     ADD_ASSIGN = 272,
     SUB_ASSIGN = 273,
     LEFT_ASSIGN = 274,
     RIGHT_ASSIGN = 275,
     AND_ASSIGN = 276,
     XOR_ASSIGN = 277,
     OR_ASSIGN = 278,
     INT = 279,
     FLOAT = 280,
     DOUBLE = 281,
     VOID = 282,
     TYPEDEF = 283,
     INLINE = 284,
     BOOL = 285,
     IF = 286,
     ELSE = 287,
     RETURN = 288,
     SIGN = 289,
     ABS = 290,
     SQRT = 291,
     COMPARE = 292,
     ENUM = 293,
     EXTERN = 294,
     EXACT = 295,
     GROUP = 296,
     DEGREE = 297,
     INT_CONSTANT = 298,
     FLOAT_CONSTANT = 299,
     IDENTIFIER = 300,
     STRING_LITERAL = 301,
     USER_TYPE = 302
   };
#endif
/* Tokens.  */
#define PTR_OP 258
#define INC_OP 259
#define DEC_OP 260
#define LEFT_OP 261
#define RIGHT_OP 262
#define LE_OP 263
#define GE_OP 264
#define EQ_OP 265
#define NE_OP 266
#define AND_OP 267
#define OR_OP 268
#define MUL_ASSIGN 269
#define DIV_ASSIGN 270
#define MOD_ASSIGN 271
#define ADD_ASSIGN 272
#define SUB_ASSIGN 273
#define LEFT_ASSIGN 274
#define RIGHT_ASSIGN 275
#define AND_ASSIGN 276
#define XOR_ASSIGN 277
#define OR_ASSIGN 278
#define INT 279
#define FLOAT 280
#define DOUBLE 281
#define VOID 282
#define TYPEDEF 283
#define INLINE 284
#define BOOL 285
#define IF 286
#define ELSE 287
#define RETURN 288
#define SIGN 289
#define ABS 290
#define SQRT 291
#define COMPARE 292
#define ENUM 293
#define EXTERN 294
#define EXACT 295
#define GROUP 296
#define DEGREE 297
#define INT_CONSTANT 298
#define FLOAT_CONSTANT 299
#define IDENTIFIER 300
#define STRING_LITERAL 301
#define USER_TYPE 302




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 56 "parse.y"
{
  AST::Expression* exp;
  AST::Statement* stm;
  AST::CompoundStatement* compound_stm;
  AST::ExpressionStatement* exp_stm;
  AST::StatementList* statements;
  AST::VariableDeclaration *var_decl;
  AST::StatementList* stm_list;


  //AST::SwitchStatement* switch_stm;
  //AST::WhileLoop* while_loop;
  //AST::DoLoop* do_loop;
  //AST::ForLoop* for_loop;
  //AST::LabelContext* label_context;

  Declarator *declarator;
  FunctionDeclarator::ParameterDeclList *param_decl_list;
  FunctionDeclarator::ParameterDecl *param_decl;
  std::list< Declarator* > *decl_list;
  AST::ExpressionList *exp_list;

  AST::BinaryExpression::Kind binexp_kind;
  AST::UnaryExpression::Kind unexp_kind;
  AST::AssignmentExpression::Kind ass_exp_kind;
  AST::FunctionDefinition *fun_def;
  //AST::SwitchStatement::SwitchCase *switch_case;
  //AST::SwitchStatement::CaseContainer *case_container;
  FunctionType *fun_type;
  //StructType *struct_type;

  Type* type;
  int int_const;
  double float_const;
  char* string_const;

  Group_rep* single_group;
  std::list< Group_rep* > *group_list;
}
/* Line 1489 of yacc.c.  */
#line 183 "parse.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

