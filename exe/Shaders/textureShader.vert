varying vec4 position;  // position of the vertex (and fragment) in world space
varying vec3 varyingNormalDirection;  // surface normal vector in world space
varying vec4 texCoords; // the texture coordinates
 
void main()
{
  position = gl_ModelViewProjectionMatrix * gl_Vertex;
  varyingNormalDirection = normalize(gl_NormalMatrix * gl_Normal); 
  texCoords = gl_MultiTexCoord0;
  gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}