var current = 0;

function next(){ // forward one image
  if(document.formname.slide[current+1]){
    document.images.show.src = document.formname.slide[current+1].value;
    document.formname.slide.selectedIndex = ++current;}
  else{first();}}

function previous(){ // back on image
  if((current-1) >= 0){
    document.images.show.src = document.formname.slide[current-1].value;
    document.formname.slide.selectedIndex= --current;}
  else{last();}}

function first(){ // jump to first image
  current=0;
  document.images.show.src = document.formname.slide[0].value;
  document.formname.slide.selectedIndex=0;}

function last(){ // this is jump to last image
  current=(document.formname.slide.length-1);
  document.images.show.src = document.formname.slide[current].value;
  document.formname.slide.selectedIndex=current;}

function ap(text){ // this controls the auto-play and/or auto-stop
  document.formname.slidebutton.value=(text == "Stop") ? "Start" : "Stop";
  rotate();}

function change(){ // this is for the pulldown menu
  current=document.formname.slide.selectedIndex;
  document.images.show.src = document.formname.slide[current].value;}

function rotate() {
  if (document.formname.slidebutton.value == "Stop") {
    current = (current == document.formname.slide.length-1) ? 0 : current+1;
    document.images.show.src = document.formname.slide[current].value;
    document.formname.slide.selectedIndex = current;
    window.setTimeout("rotate()", 2000);}}
