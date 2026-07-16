/* ============ STATISTICAL CORE (identical to what will be embedded in index.html) ============ */

function gammln(xx){
  const cof=[76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5];
  let x=xx, y=xx;
  let tmp=x+5.5; tmp-=(x+0.5)*Math.log(tmp);
  let ser=1.000000000190015;
  for(let j=0;j<6;j++){ y+=1; ser+=cof[j]/y; }
  return -tmp+Math.log(2.5066282746310005*ser/x);
}
function gser(a,x){
  const ITMAX=200, EPS=3e-9;
  const gln=gammln(a);
  if(x<=0) return {gamser:0, gln};
  let ap=a, sum=1/a, del=sum;
  for(let n=1;n<=ITMAX;n++){ ap+=1; del*=x/ap; sum+=del; if(Math.abs(del)<Math.abs(sum)*EPS) break; }
  return {gamser: sum*Math.exp(-x+a*Math.log(x)-gln), gln};
}
function gcf(a,x){
  const ITMAX=200, EPS=3e-9, FPMIN=1e-300;
  const gln=gammln(a);
  let b=x+1-a, c=1/FPMIN, d=1/b, h=d;
  for(let i=1;i<=ITMAX;i++){
    const an=-i*(i-a); b+=2;
    d=an*d+b; if(Math.abs(d)<FPMIN) d=FPMIN;
    c=b+an/c; if(Math.abs(c)<FPMIN) c=FPMIN;
    d=1/d; const del=d*c; h*=del;
    if(Math.abs(del-1)<EPS) break;
  }
  return {gammcf: Math.exp(-x+a*Math.log(x)-gln)*h, gln};
}
function gammq(a,x){
  if(x<0||a<=0) return NaN;
  if(x<a+1) return 1-gser(a,x).gamser;
  return gcf(a,x).gammcf;
}
function chiSqPValue(chi2, df){ return gammq(df/2, chi2/2); }

function chiSquareTest(table){
  const rows=table.length, cols=table[0].length;
  const rowSums=table.map(r=>r.reduce((a,b)=>a+b,0));
  const colSums=table[0].map((_,c)=>table.reduce((a,r)=>a+r[c],0));
  const N=rowSums.reduce((a,b)=>a+b,0);
  const expected=table.map((row,i)=>row.map((_,j)=>rowSums[i]*colSums[j]/N));
  let chi2=0;
  for(let i=0;i<rows;i++) for(let j=0;j<cols;j++){ const e=expected[i][j]; if(e>0) chi2+=(table[i][j]-e)**2/e; }
  const df=(rows-1)*(cols-1);
  return {chi2, df, p: chiSqPValue(chi2, df), expected};
}

function rankArray(values){
  const idx=values.map((v,i)=>i).sort((a,b)=>values[a]-values[b]);
  const ranks=new Array(values.length);
  let i=0;
  while(i<idx.length){
    let j=i;
    while(j+1<idx.length && values[idx[j+1]]===values[idx[i]]) j++;
    const avgRank=(i+1+j+1)/2;
    for(let k=i;k<=j;k++) ranks[idx[k]]=avgRank;
    i=j+1;
  }
  return ranks;
}
function kruskalWallis(groups){
  const all=groups.flat();
  const N=all.length;
  const ranks=rankArray(all);
  let offset=0; const Rsum=[];
  groups.forEach(g=>{ let s=0; for(let i=0;i<g.length;i++) s+=ranks[offset+i]; Rsum.push(s); offset+=g.length; });
  let H=(12/(N*(N+1)))*groups.reduce((acc,g,i)=>acc+(Rsum[i]**2)/g.length,0) - 3*(N+1);
  const counts={};
  all.forEach(v=>counts[v]=(counts[v]||0)+1);
  const tieSum=Object.values(counts).reduce((acc,t)=>acc+(t**3-t),0);
  const C=1-tieSum/(N**3-N);
  if(C>0) H=H/C;
  const df=groups.length-1;
  return {H, df, p: chiSqPValue(H, df)};
}

function linearRegressionStats(x,y){
  const n=x.length;
  const mx=x.reduce((a,b)=>a+b,0)/n, my=y.reduce((a,b)=>a+b,0)/n;
  let num=0, den=0;
  for(let i=0;i<n;i++){ num+=(x[i]-mx)*(y[i]-my); den+=(x[i]-mx)**2; }
  const slope=den===0?0:num/den;
  const intercept=my-slope*mx;
  let ssRes=0, ssTot=0;
  for(let i=0;i<n;i++){ const pred=intercept+slope*x[i]; ssRes+=(y[i]-pred)**2; ssTot+=(y[i]-my)**2; }
  const r2=ssTot===0?1:1-ssRes/ssTot;
  return {slope, intercept, r2};
}

function medianOf(arr){ const s=[...arr].sort((a,b)=>a-b); const n=s.length; const mid=Math.floor(n/2); return n%2===0?(s[mid-1]+s[mid])/2:s[mid]; }
function quartiles(values){
  const sorted=[...values].sort((a,b)=>a-b);
  const n=sorted.length, mid=Math.floor(n/2);
  const lower = sorted.slice(0,mid);
  const upper = n%2===0 ? sorted.slice(mid) : sorted.slice(mid+1);
  const q1=medianOf(lower), q3=medianOf(upper);
  return {q1, q3, iqr:q3-q1};
}

/* ============ VALIDATION SUITE ============ */
function approx(a,b,tol){ return Math.abs(a-b)<=tol; }
let pass=0, fail=0;
function check(name, cond, detail){
  if(cond){ pass++; console.log('PASS  -', name); }
  else { fail++; console.log('FAIL  -', name, '|', detail); }
}

console.log('=== 1. IQR / QUARTIL ===');
{
  const r = quartiles([1,2,3,4,5,6,7,8]);
  check('Q1 = 2.5 (median bawah [1,2,3,4])', approx(r.q1,2.5,1e-9), r.q1);
  check('Q3 = 6.5 (median atas [5,6,7,8])', approx(r.q3,6.5,1e-9), r.q3);
  check('IQR = 4', approx(r.iqr,4,1e-9), r.iqr);
}

console.log('\n=== 2. REGRESI LINEAR (dihitung manual: slope=0.9, intercept=1.3, R2=0.81) ===');
{
  const r = linearRegressionStats([1,2,3,4,5],[2,3,5,4,6]);
  check('slope = 0.9', approx(r.slope,0.9,1e-9), r.slope);
  check('intercept = 1.3', approx(r.intercept,1.3,1e-9), r.intercept);
  check('R2 = 0.81', approx(r.r2,0.81,1e-9), r.r2);

  const perfect = linearRegressionStats([1,2,3,4,5],[2,4,6,8,10]);
  check('kasus sempurna: slope=2, intercept=0, R2=1', approx(perfect.slope,2,1e-9)&&approx(perfect.intercept,0,1e-9)&&approx(perfect.r2,1,1e-9), JSON.stringify(perfect));
}

console.log('\n=== 3. DISTRIBUSI CHI-SQUARE (cek terhadap nilai kritis baku yang sudah dikenal luas, p seharusnya ~0.05) ===');
{
  const knownCritical = [[1,3.841],[2,5.991],[3,7.815],[4,9.488],[5,11.070],[10,18.307]];
  knownCritical.forEach(([df,x])=>{
    const p = chiSqPValue(x, df);
    check(`df=${df}, x=${x} -> p~0.05`, approx(p,0.05,0.001), p);
  });
}

console.log('\n=== 4. CHI-SQUARE INDEPENDENCE TEST (tabel 2x2 dihitung manual: chi2=0.7937, df=1) ===');
{
  const r = chiSquareTest([[10,20],[30,40]]);
  check('chi2 = 0.7937', approx(r.chi2,0.7937,0.001), r.chi2);
  check('df = 1', r.df===1, r.df);
  check('expected[0][0] = 12', approx(r.expected[0][0],12,1e-9), r.expected[0][0]);
  console.log('    p-value =', r.p, '(bandingkan manual ~0.373 via 2*(1-Phi(0.891)))');
}

console.log('\n=== 5. KRUSKAL-WALLIS (3 grup tanpa ties, solusi closed-form df=2: H=7.2, p=e^-3.6=0.0273237...) ===');
{
  const r = kruskalWallis([[1,2,3],[4,5,6],[7,8,9]]);
  check('H = 7.2', approx(r.H,7.2,1e-9), r.H);
  check('df = 2', r.df===2, r.df);
  check('p = exp(-3.6) = 0.0273237...', approx(r.p, Math.exp(-3.6), 1e-6), r.p);
}

console.log('\n=== 6. KRUSKAL-WALLIS dengan TIES (cek koreksi tie tidak merusak hasil - grup identik harus H=0) ===');
{
  const r = kruskalWallis([[5,5,5],[5,5,5],[5,5,5]]);
  check('semua nilai identik -> H = 0', approx(r.H,0,1e-6), r.H);
}

console.log(`\n=== HASIL: ${pass} PASS, ${fail} FAIL ===`);
process.exit(fail>0?1:0);
