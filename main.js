function logit(text) {
  if(window.console)
    console.log(text);
}

$(document).ready(main);
var cqNoise = new CqNoise();

function main() {
  var canvas = $('#canvas')[0];
  context = canvas.getContext('2d');
  context.fillStyle = 'black';
  context.fillRect(0, 0, canvas.width, canvas.height);
  var imageData = context.getImageData(0, 0, canvas.width, canvas.height);
  var x, y;
  var scale = 0.013;
  for(x = 0; x < canvas.width; x++) {
    for(y = 0; y < canvas.height; y++) {
      imageData.data[(y * canvas.width + x) * 4 + 0] = 0;
      imageData.data[(y * canvas.width + x) * 4 + 1] = (cqNoise.noise2(x * scale, y * scale) * 0.5 + 0.5) * 255;
      imageData.data[(y * canvas.width + x) * 4 + 2] = 0;
      imageData.data[(y * canvas.width + x) * 4 + 3] = 255;
    }
  }
  context.putImageData(imageData, 0, 0);
  logit('painted');
}
