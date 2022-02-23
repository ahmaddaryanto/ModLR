* Encoding: UTF-8.
/*Moderated Linear Regression - ModLR. 
/*macro created by Ahmad Daryanto (a.daryanto@lancaster.ac.uk)*/. 

set printback=off mprint=off errors=off .
preserve .
* Encoding: UTF-8.
SET MXLOOPS = 10000001.
set printback off.
DEFINE MLR (iv = !charend('/')
/covs = !charend('/') 
/mod = !charend('/')
/dv = !charend('/')
/tVar = !charend('/') !default(0) 
/tLine = !charend('/') !default(0)
/tiv = !charend('/') !default(0)
/tcat = !charend('/') !default(0)
/hci = !charend('/') !default(0)
/ramsey = !charend('/') !default(0)
/quadratic = !charend('/') !default(0)
/multim = !charend('/') !default(0))    
MATRIX.
get mat/variables=!covs !iv !mod !dv /names=nms /MISSING=OMIT.
compute n=nrow(mat).
compute k=ncol(mat).
compute ones=make(n,1,1).
*dv in original metrix.
compute Y=mat(:,k).
compute notft=1.
*=======================.
*         Descriptive statistics.
*======================.
print /title="ModLR (v250619) written by Ahmad Daryanto".
print /title="https://sites.google.com/site/ahmaddaryanto/scripts/ModLR".
compute temp1=t(nms(:,k)).
compute temp2=t(nms(:,(k-2))).
compute temp3=t(nms(:,(k-1))).
compute cnms={"DV ==>","IV ==>", "Mod ==>"}.
compute temp=t({temp1,temp2,temp3}).
print temp/title="Key variables"
/rnames=cnms
/ format=A10.
print n/title " Sample size ".
compute checkint=1.
loop #i=1 to nrow(mat(:,k-1)).
compute modtes=mod(mat(:,k-1)(#i),1).
do if modtes NE 0.
compute checkint=0.
end if.
end loop.   
do if (!tcat=1) AND checkint=1.
*compute subgroup sample size (group a and b). 
compute n_grp={0,0}.
loop #i=1 to n.
do if mat(#i,(k-1))=0.  
compute n_grp(1)=n_grp(1)+1.     
else if mat(#i,(k-1))=1.
compute n_grp(2)=n_grp(2)+1.   
end if.
end loop.
compute cnms={"group 0=", "group 1="}.
print t(n_grp)/title= "Moderator is a (0,1) variable, within group sample sizes:"
/rnames=cnms.
else if (!tcat=1) AND checkint=0.
print/title='WARNING:'.
print/title='Moderator is specified as a categorical variable but not all values are '+
'integer'.
print/title='------------------------------------'.
end if.
*compute temp=t(nms(:,1:(k-3))).
*print temp/title="Control Variable(s)"
*/format=A8.
*mean, sd calculation.
compute Xbar=(1/n)*ones*t(ones)*mat.
compute D=mat-Xbar. 
compute CM = 1/n*sscp(D).
compute V = diag(CM).
compute sd=diag((mdiag(sqrt(V)))).
compute S = inv(mdiag(sqrt(V))).
compute mean=t(Xbar(1,:)).
* computing corr matrix.
compute R = S*CM*S.
compute cnms={"Mean","SD", "Corr"}.
compute nmvars = t(nms(1,1:k)).
compute nmvars = {nmvars}.
compute desc={mean,sd,R}.
print desc/title "Mean, SD, Correlation"/rnames=nmvars/cnames=cnms/format f9.3.
*--check determinant.
compute X={ones,mat(:,1:(k-1))}.
do if det(sscp(X)) <= 0.
print /title="========================================".
print /title="Non-singular matrix. Matrix is invertible. Program Stops.".
print /title="========================================".
else if det(sscp(X)) > 0.
*====== Prepare the OLS matrix, eg. 0=original, 1=centered ======.
do if (!tVar=1).                          
compute X={mat(:,k-2), mat(:,k-1)}. 
compute Xbar=(1/n)*ones*t(ones)*X.
compute D=X-Xbar.                    
compute inter=D(:,1) &* D(:,2).
compute mtrick={ones,mat}.
compute mshadow=mtrick(:,1:(k-2)).
compute X={mshadow,D,inter}.
*save a constant and covariates-are not centered.
compute Xcov=X(:,1:(k-2)).
else if (!tVar=0).
* the matrix without Y.           
compute inter=mat(:,(k-2))&*mat(:,(k-1)).
compute X={ones,mat(:,1:(k-1)),inter}.
*save a constant and covariates-are not centered.
compute Xcov=X(:,1:(k-2)).
end if. 
*===============================================.
*            OLS Regression .
*===============================================.
compute b=(inv(sscp(X)))*t(X)*Y.  
*===computing standard error of b, t value and p-value of 
OLS centered==.
*****se(b).Verbeek p.19.
compute j=ncol(X).
compute e=Y-X*b.
compute sserr=sscp(e).
compute mse=(1/(n-j))*sserr.
*--robust standard error.
do if (!hci=0).
compute vb=mse*inv(sscp(X)).
else if (!hci=1).
compute esq=e&*e.
compute mesq=mdiag(esq).
compute  vb=inv(sscp(X))*t(X)*mdiag(esq)*X*inv(sscp(X)).
end if.  
compute sb=sqrt(diag(vb)).
compute tb=b/sb.
compute dff=n-j.
compute F=tb&*tb.
compute pF=1-fcdf(F,1,dff).
*--95% CI--.
compute LB=b-1.96*sb.
compute UB=b+1.96*sb.
*===Preparing input ANOVA table.
*computing mean square regression.
compute meanY=ones*mean(k).
compute e_reg=X*b-meanY.
compute ssreg=csum(sscp(e_reg)).
compute sumsq=T({ssreg,sserr}).
compute dfa=T({k,n-k-1}).
compute mse_a=sumsq/dfa.
compute Fval=(ssreg/(k))/(sserr/(n-k-1)).
Compute pF_a=1-fcdf(Fval,k,n-k-1).
compute F_a=T({Fval,-999}).
compute pFa=T({pF_a,-999}).
*----Akaike----.
compute aic_res=n*ln(sserr/n)+(2*j)+(n*(1+ln (2*3.14159))).
*===============================================.
*===computing R-sq for regression with an interaction term.
compute M = (1/n)*ones*t(ones)*mat.
compute My ={M(:,k)}. 
compute dY = y- My.
compute sum_dY=(1/(n-1))*csum(sscp(dY)).
compute sum_e2=(1/(n-1))*csum(sscp(e)).
compute Rsq_in=1-sum_e2/sum_dY.
*===============================================.
*            OLS Regression of without interaction.
*===============================================.
compute Xu={ones,mat(:,1:(k-1))}.
compute ju=ncol(Xu).
compute bu=(inv(sscp(Xu)))*t(Xu)*Y.
*SSer=sum of squares error.
compute eu=Y-Xu*bu.
compute sseru=sscp(eu).
compute mseu=(1/(n-1))*sseru.
compute Rsq_out=1-mseu/sum_dY.
*--robust standard error.
do if (!hci=0).
compute vbu=mse*inv(sscp(X)).
else if (!hci=1).
compute esq=e&*e.
compute mesq=mdiag(esq).
compute vbu=inv(sscp(X))*t(X)*mdiag(esq)*X*inv(sscp(X)).
end if.  
*Akaike----.
compute aic_main=n*ln(sseru/n)+(2*ju)+(n*(1+ln (2*3.14159))).
*================================================.
*===computing F test for significance of the interaction term.
*Verbeek p.59.
compute w=ncol(X)-ncol(Xu).
compute dfnum=w.
compute dfdenom=n-k-1.
compute Fc=((Rsq_in-Rsq_out)/(dfnum))/((1-Rsq_in)/(dfdenom)).
compute pFc=1-fcdf(Fc,dfnum,dfdenom).
compute diff=Rsq_in-Rsq_out.
compute Rsq=T({Rsq_out, Rsq_in}).
compute change={diff, dfnum,dfdenom, Fc, pFc}.
*computing effect size (Cohen 1988, Aiken & West, p.157).
compute Effsz=diff/(1-Rsq_in).
*================================================.
*==== RAMSEY'S RESET TEST ======*.
do if (!ramsey=1).
compute yhat=X*b.
compute yhat2=yhat &* yhat.                                             
compute Xr={X,yhat2}.
compute br=(inv(sscp(Xr)))*t(Xr)*Y. 
compute jr=ncol(Xr).
compute e=Y-Xr*br.
compute sser=sscp(e).
compute mse=(1/(n-jr))*sser.
*--robust standard error.
do if (!hci=0).
compute vb=mse*inv(sscp(Xr)).
else if (!hci=1).
compute esq=e&*e.
compute mesq=mdiag(esq).
compute  
vb=inv(sscp(Xr))*t(Xr)*mdiag(esq)*Xr*inv(sscp(Xr)).
end if.  
compute sbr=sqrt(diag(vb)).
compute tbr=br/sbr.
compute dffr=n-jr.
compute Fr=tbr&*tbr.
compute pFr=1-fcdf(Fr,1,dffr).                              
end if. 
*=== end of RESET===.
*====F test for joint restrictions on the quadratic terms======.
do if (!quadratic=1).
*---compute sum of squares for unresticted model, i.e., with quadratic 
terms.
compute x1=mat(:,k-2).
compute z=mat(:,k-1).
compute xsq=x1 &* x1.
compute zsq=z &* z.
compute Xur={X,xsq,zsq}.
compute jur=ncol(Xur).
do if det(sscp(Xu))>0.    
compute NotFt=0.                                                    
compute bur=(inv(sscp(Xur)))*t(Xur)*Y.
compute e=Y-Xur*bur.
compute sser_ur=sscp(e).
compute q=2.      
compute dfq=n-jur. 
compute Fquad=((sserr-sser_ur)/q)/(sser_ur/dfq).                    
compute pFquad=1-fcdf(Fquad,q,dfq).  
*----Akaike----.
compute aic_quad=n*ln(sser_ur/n)+(2*jur)+(n*(1+ln (2*3.14159))).    
RELEASE x1.
RELEASE z.
compute mse=(1/(n-jur))*sser_ur.                                 
*--robust standard error.
do if (!hci=0).
compute vbur=mse*inv(sscp(Xur)).
else if (!hci=1).
compute esq=e&*e.
compute mesq=mdiag(esq).
compute 
vbur=inv(sscp(Xur))*t(Xur)*mdiag(esq)*Xur*inv(sscp(Xur)).
end if.  
compute sbur=sqrt(diag(vbur)).
compute tbur=bur/sbur.                                          
compute F=tbur&*tbur.
compute pFur=1-fcdf(F,1,dfq).   
else if det(sscp(Xu))<=0.  
compute NotFt=1.
end if.                                                             
end if.
*======= end of for joint restrictions on the quadratic terms=====.
*=============================================*.
*                         Simple Effect Analysis  *.
*=============================================*.
compute btrick=t(b((j-2):j)).
compute bsef=t({b(1),btrick}).
do if (!tVar=1).
compute ivlo=-sd(k-2).
compute ivhi=sd(k-2).
do if (!tcat=1).
compute modlo=0.
compute modhi=1.
else if (!tcat=0).
compute modlo=-sd(k-1).
compute modhi=sd(k-1).
end if. 
compute modmean=0.
compute ivmean=0.
else if (!tVar=0).
compute ivlo=mean(k-2)-sd(k-2).
compute ivhi=mean(k-2)+sd(k-2).
do if (!tcat=1).
compute modlo=0.
compute modhi=1.
else if (!tcat=0).
compute modlo=mean(k-1)-sd(k-1).
compute modhi=mean(k-1)+sd(k-1).
end if.
compute modmean=mean(k-1).
compute ivmean=mean(k-2).
end if. 
*--the plot option--- .
do if (!tLine=0).
compute modr={modlo;modhi;modlo;modhi}.
compute iv_={ivlo;ivlo;ivhi;ivhi}.
else if (!tLine=1).
compute modr=
{modlo;modmean;modhi;modlo;modmean;modhi;modlo;modmean;modhi}.
compute iv_={ivlo;ivlo;ivlo;ivmean;ivmean;ivmean;ivhi;ivhi;ivhi}.
end if. 
compute dv_=bsef(1)+bsef(2)*iv_+bsef(3)*modr+bsef(4)*iv_&*modr.
compute dataplot={iv_,modr,dv_}.
save dataplot/outfile=* /var  IV M DV.
*--------------- computing beta and se of the slopes -----------------.
do if (!tcat=0).
compute bl={999,999,999}.
compute sl={999,9999,999}.
compute tl={999,9999,999}.
compute pl={999,9999,999}.
compute lbl={999,9999,999}.
compute ubl={999,9999,999}.
else if (!tcat=1).
compute bl={999,999}.
compute sl={999,9999}.
compute tl={999,9999}.
compute pl={999,9999}.
compute lbl={999,9999}.
compute ubl={999,9999}.
end if.
*--- op contains mean-1SD,Mean, and mean+1SD.
do if (!tcat=1).
compute op={0,1}.
else if (!tcat=0) and (!tVar=1).                                         
compute op={-sd(k-1),0,sd(k-1)}.      
else if (!tcat=0) and (!tVar=0).                                         
compute op={mean(k-1)-sd(k-1),mean(k-1),mean(k-1)+sd(k-1)}.    
end if.
do if (!tVar=0).
compute iv=mat(:,k-2).
compute z=mat(:,k-1).
else if (!tVar=1).
compute iv=D(:,1).       
compute z=D(:,2).                    
end if.
/*preparing X matrix for simple effect.
do  if ncol(Xcov)=1.
compute Xsimple={ones,iv}.
else if ncol(Xcov)>1.
compute Xsimple={Xcov,iv}.
end if.
do if (!tcat=0).
   compute numloop=3.
else if (!tcat=1).
  compute numloop=2.
end if.
loop #i=1 to numloop. 
compute at_z=z-op(#i).  
compute inter=iv &* at_z.                                 
compute X={Xsimple,at_z,inter}.
compute bs=(inv(sscp(X)))*t(X)*Y. 
compute j=ncol(X).
compute e=Y-X*bs.
compute e2=sscp(e).
compute sser=csum(e2).
compute mse=(1/(n-j))*sser.
*--robust standard error.
do if (!hci=0).
compute vs=mse*inv(sscp(X)).
else if (!hci=1).
compute esq=e&*e.
compute mesq=mdiag(esq).
compute  vs=inv(sscp(X))*t(X)*mdiag(esq)*X*inv(sscp(X)).
end if.  
compute ss=sqrt(diag(vs)).
compute ts=bs/ss.
compute dff=n-j.
compute Fs=ts&*ts.
compute ps=1-fcdf(Fs,1,dff).
*--95% CI--.
compute lbs=bs-1.96*ss.
compute ubs=bs+1.96*ss.
compute gotit=nrow(bs)-2.
compute bl(#i)=bs(gotit).
compute sl(#i)=ss(gotit).
compute tl(#i)=ts(gotit).
compute pl(#i)=ps(gotit).
compute lbl(#i)=lbs(gotit).
compute ubl(#i)=ubs(gotit).
end loop.                         
*this part was adapted from O'Connor (1998)
* df & t values; from Darlington p 516 & Howell 87 p 586; p = 05 two-tailed 
.
compute dft={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,    
30,32,34,36,38,40,43,46,49,52,56,60,65,70,75,80,85,90,95,100,110,120,130,
150,175,200,250,300,400,500,600,700,800,900,1000,1000000000;
12.706,4.303,3.182,2.776,2.571,2.447,2.365,2.306,2.262,2.228,2.201,2.179,
2.160,2.145,2.131,2.120,2.110,2.101,2.093,2.086,2.074,2.064,2.056,2.048,
2.042,2.037,2.032,2.028,2.024,2.021,2.017,2.013,2.010,2.007,2.003,2.000,
1.997,1.994,1.992,1.990,1.988,1.987,1.985,1.984,1.982,1.980,1.978,1.976,
1.974,1.972,1.969,1.968,1.966,1.965,1.964,1.963,1.963,1.963,1.962,1.962 }.
compute tabledT = 0.
loop #a = 1 to 59  .
do if (dfdenom ge dft(1,#a) and dfdenom < dft(1,#a+1) ).
compute tabledT = dft(2,#a) .
end if.
end loop if (tabledT > 0).
*=========================================================.
*==Johnson-Neymann Technique for Calculation Limits of Insignificance 
Region.
*    see Huitema p.249 and Pedhazur 97 p.469.
*    Moderator is categorical variable.
*=========================================================.
*            OLS Regression of Raw Data for Each Group.
*=========================================================.
compute checkint=1.
loop #i=1 to nrow(mat(:,k-1)).
compute modtes=mod(mat(:,k-1)(#i),1).
do if modtes NE 0.
compute checkint=0.
end if.
end loop.
print/title='------------------------------------'.
do if (!tcat=1) AND checkint=0.
print/title='WARNING:'.
print/title='Moderator is specified as a categorical but not all '+
'values are integer'.
print/title='Simple slope is calculated at M=0 and M=1'.
print/title='------------------------------------'.
end if.
do if (!tcat = 1) AND (checkint=1).     
do if pF(k+1) LT 0.0500001.
compute flag=1.
*create matrix containing data for each group.
compute matg_a=make(n_grp(1),ncol(mat),modlo).
compute matg_b=make(n_grp(2),ncol(mat),modhi).
compute i=1.
compute j=1.
loop #k=1 to nrow(mat).
do if  mat(#k,(k-1))=modlo.    
compute matg_a(i,1)=mat(#k,k).
compute matg_a(i,2)=mat(#k,1).    
compute i=i+1.
else if  mat(#k,(k-1))=modhi.    
compute matg_b(j,1)=mat(#k,k).
compute matg_b(j,2)=mat(#k,1).    
compute j=j+1.
end if.
end loop. 
*group regression.
*group a.
compute ones=make(nrow(matg_a),1,1).
compute Y=matg_a(:,1).
compute X={ones,matg_a(:,2:2)}.
compute ba=(inv(sscp(X)))*t(X)*Y.
compute e=Y-X*ba.     
compute e2=sscp(e). 
compute sres_a=csum(e2).
*group b.    
compute ones=make(nrow(matg_b),1,1).
compute Y=matg_b(:,1).
compute X={ones,matg_b(:,2:2)}.
compute bb=(inv(sscp(X)))*t(X)*Y.
compute e=Y-X*bb.     
compute e2=sscp(e). 
compute sres_b=csum(e2).
*pooled residual sum of squares.
compute ssres=sres_a+sres_b.
*mean group and sum of squared deviation.
compute x_sq={1,2}.
compute m_grp={1,2}.
compute M = 
csum(matg_a(:,2))*1/nrow(matg_a(:,2)).
compute D = matg_a(:,2)-M.
compute x_sq(1)=mssq(D).
compute m_grp(1)=M.
compute M = 
csum(matg_b(:,2))*1/nrow(matg_b(:,2)).
compute D = matg_b(:,2)-M.
RELEASE matg_a.
RELEASE matg_b.
RELEASE e.
RELEASE y.
compute x_sq(2)=mssq(D).
compute m_grp(2)=M.
*computing F(alpha).
compute dfval=n-4.
compute tabledT = 0.
loop #a = 1 to 59  .
do if (dfval GE dft(1,#a) and dfval LT 
dft(1,#a+1) ).
compute tabledT = dft(2,#a) .
end if.
end loop if (tabledT > 0).
compute F_jn=tabledT**2.
compute 
Aconst=(-F_jn/dfval)*ssres*(1/(x_sq(1))+1/(x_sq(2)))+(bb(2) -ba(2))**2.
compute 
Bconst=(F_jn/dfval)*ssres*(m_grp(1)/(x_sq(1))+m_grp(2)/(x_sq(2)))+(bb(2)-ba(2))*(bb(1)-ba(1)).
compute 
Cconst=(-F_jn/dfval)*ssres*(n/(n_grp(1)*n_grp(2))+m_grp(1)** 2/(x_sq(1))+m_grp(2) 
**2/(x_sq(2)))+(bb(1)-ba(1))**2.
compute 
LowCI=(-Bconst-sqrt(Bconst**2-Aconst*Cconst))*(1/Aconst).
compute 
UpCI=(-Bconst+sqrt(Bconst**2-Aconst*Cconst))*(1/Aconst).
do if(!tVar=1).
compute X_LB=LowCI.
compute X_UB=UpCI.
else if (!tVar=0).
compute X_LB=LowCI-mean(1).
compute X_UB=UpCI-mean(1).
end if.                
else.
compute flag=0.
end if.
else.
*print /title "* Moderator is not categorical *".
compute flag=0.
end if.
RELEASE ONES.
RELEASE X.  
RELEASE M.
RELEASE D.
*==================OLS OUTPUT===========================.
Do if (!tVar=1).
print /title '** In regression, variables of the interaction '+
'term were mean-centered **'.
else if  (!tVar=0).
print /title ' ** In regression, original variables were used**'.
end if.
compute nmvars = t(nms(1,1:(k-1))).
compute nmvars = {"constant"; nmvars; "interact"}.
compute olsOut={b,sb,tb,pF, LB,UB}.
compute cnms={"b","se", "t", "sig", "95%LB", "95%UB"}.
do if(!hci=1).
print /title="** Values for se: Heteroskedasticity-robust "+
"standard errors, variant HC1 **".
end if.
print olsOut/title "Regression, unstandardized estimates "+
""/rnames=nmvars/cnames=cnms/format F9.3.
*--=================ANOVA table=============.
print {sumsq,dfa, mse_a,F_a,pFa} /space=3
/title '------- ANOVA TABLE Regression with Interaction --------'
/clabel "SS" "df" "MS" "F" "Sig"
/rlabel "Model" "Error"
/format f10.3 .
*==================================================.
compute desc={Rsq}.
compute cnms={"R-square", "dfnum", "dfdenom"}.
compute rnms={"no inter", "inter"}.
print desc/title ="Summary"/rnames=rnms/cnames=cnms/format=F9.3.
compute cnms={"change", "dfnum", "dfdenom","F","sig"}.
compute rnms={" "}.
print change/title =" "/rnames=rnms/cnames=cnms/format=F9.3.
print {Effsz} 
/title 'Effect Size (f-square) of the Interaction Term (Cohen 1988, '+
'Aiken & West 1996, p.157) '
/format f9.3 .
*==================================================.
*=================Simple slope output ===================.
do if (!tcat=0).
compute nmvars = {"M-1SD";"M"; "M+1SD"}.
else if (!tcat=1).
compute nmvars = {"0";"1"}.
end if.
compute cnms={"at","b","se", "t", "sig", "95%LB", "95%UB"}.
compute simple={t(op),t(bl),t(sl),t(tl),t(pl), t(lbl),t(ubl)}.
print simple/title "Simple slopes, unstandardized estimates (Note: M=mean value of the moderator) "
/rnames=nmvars/cnames=cnms/format F9.3.    
*====J-N Region==========.
*===RAMSEY'S RESET TEST outputs==========.                            
do if(!ramsey=1).
compute nmvars = t(nms(1,1:(k-1))).
compute nmvars = {"constant"; nmvars; "interact";"yhat^2"}.
compute outreset={br,sbr,tbr,pFr}.
compute cnms={"b","se", "t", "sig"}.
print /title="==============================================="+
" ===================".
print outreset/title="Ramsey's RESET specification test: Auxiliary "+
"moderated regression"+
""/rnames=nmvars/cnames=cnms/format F9.3.
print /title="H0: The coefficient of yhat^2 is zero.".                      
compute dfnum=1.
compute outr={dfnum,dffr,Fr(jr), pFr(jr)}.
compute cnms={"dfnum", "dfdenom","F","sig"}.
compute rnms={" "}.
print outr/title ="Summary "/rnames=rnms/cnames=cnms/format=F9.3.
print /title=" If sig-value less than 0.05, reject the null hypothesis. ".
end if.
*======Joint F-test for quadratic terms =======================.
do if (!quadratic=1) AND (NotFt=0).  
print /title="==============================================="+
" ===================".
print /title="  Joint F-test of the significance of quadratic "+
"terms   ". 
compute nmvars = t(nms(1,1:(k-1))).
compute nmvars = {"constant"; nmvars; "interact";"X^2";"M^2"}.
compute out={bur,sbur,tbur,pFur}.
compute cnms={"b","se", "t", "sig"}.
print out/title="Regression with quadratic terms"
/rnames=nmvars/cnames=cnms/format f9.3.
compute pFquad=1-fcdf(Fquad,q,dfq).         
compute quadout={q,dfq, Fquad, pFquad}.                           
print /title="H0: The coefficients of the quadratic terms (X^2 and "+
""+
"M^2) are both zero.". 
compute cnms={"dfnum","dfdenom", "F", "sig"}.                   
print quadout/title=" Summary "/rnames=rnms/cnames=cnms/format 
f9.3. 
print /title=" If sig-value less than 0.05, reject the null "+
"hypothesis. ".
else if (!quadratic=1) AND (NotFt=1). 
print 
/title="==============================================="+
" ===================".
print /title'Joint F test was requested, but'.
print /title="Matrix is singular, thus invertible. "+
"Test can not continue. ".
end if. 
*========================End join F test=================================.
       
do if !multim=1.
 print/title="========================".
 print/title="Mult-model inference is requested".
*--Multimodel inference.
compute ones=make(n,1,1).
compute aics={0,0,0,0,0,0,0,0}.
compute jm={0,0,0,0,0,0,0,0}.
compute rss={0,0,0,0,0,0,0,0}.
loop #i=1 to 8.
do if (#i=1).
*model 1 with main effects only: x+z.
compute Xq={ones,mat(:,1:(k-1))}.
else if (#i=2).
*model 2: x+x**2.
compute X={ones,mat(:,1:(k-2))}.
compute x1=mat(:,k-2).
compute xsq=x1 &* x1.
compute Xq={ones,x1,xsq}.
else if (#i=3).
*model 3: z+z**2.
compute X={ones,mat(:,1:(k-2))}.
compute z=mat(:,k-1).
compute zsq=z &* z.
compute Xq={ones,z,zsq}.
else if (#i=4).
*model 4 with interaction: x+z+x*z.
compute inter=mat(:,(k-2))&*mat(:,(k-1)).
compute Xq={ones,mat(:,1:(k-1)),inter}.   
else if (#i=5).
*model 5: x+z+x^2+z^2.
compute X={ones,mat(:,1:(k-1))}.
compute x1=mat(:,k-2).
compute z=mat(:,k-1).
compute xsq=x1 &* x1.
compute zsq=z &* z.
compute Xq={X,xsq,zsq}.
else if (#i=6). 
*model 6: x+z+x^2+z^2+xz. 
compute inter=mat(:,(k-2))&*mat(:,(k-1)).
compute Xq={ones,mat(:,1:(k-1)),inter,xsq,zsq}.
else if (#i=7). 
*model 7: x+z+x^2+z^2+xz+z*x^2. 
compute z_xsq=z &* xsq.
compute inter=mat(:,(k-2))&*mat(:,(k-1)).
compute Xq={ones,mat(:,1:(k-1)),inter,xsq,zsq,z_xsq}.
else if (#i=8). 
*model 8: x+z+x^2+z^2+xz+x*z^2. 
compute x_zsq=x1 &* zsq.
compute inter=mat(:,(k-2))&*mat(:,(k-1)).
compute Xq={ones,mat(:,1:(k-1)),inter,xsq,zsq,x_zsq}.
end if.
compute jm(#i)=ncol(Xq)+1. 
compute b=(inv(sscp(Xq)))*t(Xq)*Y. 
compute rss(#i)=sscp(Y-Xq*b).
compute co=112/(n-8).
compute aics(#i)=n*ln(rss(#i)/n)+(2*jm(#i))+(n*(1+ln (2*3.14159)))+co.
* correction for N/K<40.
do if (n/8) <40.            
compute aics(#i)=aics(#i)+144/(n-9).
end if.
end loop.
release Xq.
compute min_aic=mmin(aics).
compute d_aic=aics-min_aic.                      
compute eps_aic={0,0,0,0,0,0,0,0}.
loop #i=1 to ncol(aics).   
compute eps_aic(#i)=2.718**(-0.5*(d_aic(#i))).                  
end loop.
compute sum_eps=rsum(eps_aic).
compute w_aic=eps_aic/sum_eps.
compute aicout={jm;rss;aics;d_aic;w_aic}.
* compute rnms={"x+x^2", "z+z^2", "x+z ", " x+z+x*z ","x+z+x^2+z^2 ", 
"x+z+x^2+z^2+xz "," x^2+z^2+x+z+xz+z*x^2","x^2+z^2+x+z+xz+x*z^2"}.
compute rnms={"1", "2","3", "4", "5", "6","7","8"}.
compute cnms={"K","RSS", "AIC", "D","Weight"}.
print /title=" Investigating Moderated Regression Models:". 
print /title="Model selection with Akaike Information Criterion (AIC)".
*print /title=" written by Ahmad Daryanto".
*print /title="sites.google.com/site/ahmaddaryanto/ ".
print /title="==================================== ".
print /title=" Models considered without writing down the model "+
"parameters and error terms:".
print /title="model 1: x+z".
print /title="model 2: x+x^2".
print /title="model 3: z+z^2".            
print /title="model 4: x+z+xz".
print /title="model 5: x^2+z^2+x+z".
print /title="model 6: x^2+z^2+x+z+xz".
print /title="model 7: x^2+z^2+x+z+xz+z*x^2".
print/title="model 8: x^2+z^2+x+z+xz+x*z^2".
print t(aicout)/title=" Summary"
/rnames=rnms
/cnames=cnms/format f9.3.
* correction for N/K<40.
do if (n/8) <40.            
print /title="AIC correction was applied because N/K is less than 40".
end if.
print/title="--------------------------".
print /title="Recommendation (Burnham & Anderson, 2003):".
print/title="D≤2     model has a substantial support.".
prin/title="4≤D≤7 model has a low support.".
print/title="D>10   model has no support.".
end if.
end if.
END MATRIX.
graph /line = mean(dv) by iv by m.
!ENDDEFINE.
restore.
